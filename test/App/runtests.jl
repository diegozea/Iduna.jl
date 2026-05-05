@testset "App" begin
    @testset "minimal CLI" begin
        kwargs = Iduna.parse_app_args(["P20963", "--mmseqs-db", "db"])
        @test kwargs[:id] == "P20963"
        @test kwargs[:mmseqs_db] == "db"
        @test kwargs[:overwrite] === false
        @test kwargs[:allow_specieslist_timeout_fallback] === true
        @test !haskey(kwargs, :pid_thresholds)
        @test !haskey(kwargs, :transcript_query_retries)
    end

    @testset "option mapping" begin
        kwargs = Iduna.parse_app_args([
            "P20963",
            "--mmseqs-db", "db",
            "--workdir", "work",
            "--output-dir", "out",
            "--overwrite",
            "--no-specieslist-timeout-fallback",
            "--threads", "4",
            "--transcript-query-retries", "3"
        ])
        @test kwargs[:workdir] == "work"
        @test kwargs[:output_dir] == "out"
        @test kwargs[:overwrite] === true
        @test kwargs[:allow_specieslist_timeout_fallback] === false
        @test kwargs[:threads] == 4
        @test kwargs[:transcript_query_retries] == 3
    end

    @testset "custom value parsing" begin
        kwargs = Iduna.parse_app_args([
            "P20963",
            "--mmseqs-db", "db",
            "--pid-thresholds", "10,20,30",
            "--transcript-query-timeout-seconds", "none",
            "--transcript-query-timeout-max-seconds", "240",
            "--thoraxe-timeout-seconds", "3600"
        ])
        @test kwargs[:pid_thresholds] == [10.0, 20.0, 30.0]
        @test kwargs[:transcript_query_timeout_seconds] === nothing
        @test kwargs[:transcript_query_timeout_max_seconds] == 240.0
        @test kwargs[:thoraxe_timeout_seconds] == 3600.0
    end

    @testset "required arguments" begin
        @test_throws Iduna.ArgParse.ArgParseError Iduna.parse_app_args(["P20963"])
        @test_throws Iduna.ArgParse.ArgParseError Iduna.parse_app_args([
            "--mmseqs-db", "db"])
    end
end
