# Julia's retry helper retries exceptions, so transient statuses are wrapped.
_is_transient_http_status(status::Integer)::Bool = status in _TRANSIENT_HTTP_STATUSES
_is_retryable_http_exception(_state, err)::Bool = err isa _RetryableHTTPStatus

_http_get_request(url::AbstractString; kwargs...) = HTTP.request("GET", url; kwargs...)

function _http_get_with_retries(url::AbstractString,
        headers;
        retries::Integer = 4,
        sleep_seconds::Real = 1.5,
        max_delay::Real = 30.0,
        http_get::Function = _http_get_request)
    attempts = max(Int(retries), 1)
    try
        # Return the final HTTP response so callers decide what status means success.
        return retry(;
            delays = Base.ExponentialBackOff(;
                n = attempts - 1,
                first_delay = Float64(sleep_seconds),
                max_delay = Float64(max_delay),
                factor = 2.0,
                jitter = 0.0),
            check = _is_retryable_http_exception) do
            resp = http_get(url; headers = headers, retry = false, status_exception = false)
            _is_transient_http_status(resp.status) && throw(_RetryableHTTPStatus(resp))
            return resp
        end()
    catch err
        err isa _RetryableHTTPStatus && return err.response
        rethrow()
    end
end
