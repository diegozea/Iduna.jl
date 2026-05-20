@testset "App" begin
    @testset "minimal CLI" begin
        kwargs = Iduna._parse_app_args(["P20963", "--mmseqs-db", "db"])
        @test kwargs[:id] == "P20963"
        @test kwargs[:mmseqs_db] == "db"
        @test kwargs[:overwrite] === false
        @test kwargs[:allow_specieslist_timeout_fallback] === true
        @test kwargs[:specieslist_filter] === true
        @test kwargs[:biomart_datasets_filter] === true
        @test kwargs[:centroids] === false
        @test kwargs[:no_expansion] === false
        @test !haskey(kwargs, :pid_thresholds)
        @test !haskey(kwargs, :transcript_query_retries)
    end

    @testset "no expansion CLI" begin
        kwargs = Iduna._parse_app_args(["P20963", "--no-expansion"])
        @test kwargs[:id] == "P20963"
        @test kwargs[:no_expansion] === true
        @test !haskey(kwargs, :mmseqs_db)
    end

    @testset "option mapping" begin
        kwargs = Iduna._parse_app_args([
            "P20963",
            "--mmseqs-db", "db",
            "--workdir", "work",
            "--output-dir", "out",
            "--overwrite",
            "--no-specieslist-timeout-fallback",
            "--orthology", "1:n",
            "--no-specieslist-filter",
            "--no-biomart-datasets-filter",
            "--centroids",
            "--threads", "4",
            "--transcript-query-retries", "3"
        ])
        @test kwargs[:workdir] == "work"
        @test kwargs[:output_dir] == "out"
        @test kwargs[:overwrite] === true
        @test kwargs[:allow_specieslist_timeout_fallback] === false
        @test kwargs[:orthology] == "1:n"
        @test kwargs[:specieslist_filter] === false
        @test kwargs[:biomart_datasets_filter] === false
        @test kwargs[:centroids] === true
        @test kwargs[:threads] == 4
        @test kwargs[:transcript_query_retries] == 3
    end

    @testset "custom value parsing" begin
        kwargs = Iduna._parse_app_args([
            "P20963",
            "--mmseqs-db", "db",
            "--pid-thresholds", "10,20,30",
            "--transcript-query-timeout-seconds", "none",
            "--transcript-query-timeout-max-seconds", "240",
            "--thoraxe-timeout-seconds", "3600",
            "--pid-sample-count", "12",
            "--pid-sample-fraction", "0.65",
            "--pid-sample-seed", "42"
        ])
        @test kwargs[:pid_thresholds] == [10.0, 20.0, 30.0]
        @test kwargs[:transcript_query_timeout_seconds] === nothing
        @test kwargs[:transcript_query_timeout_max_seconds] == 240.0
        @test kwargs[:thoraxe_timeout_seconds] == 3600.0
        @test kwargs[:pid_sample_count] == 12
        @test kwargs[:pid_sample_fraction] == 0.65
        @test kwargs[:pid_sample_seed] == 42
    end

    @testset "required arguments" begin
        @test_throws Iduna.ArgParse.ArgParseError Iduna._parse_app_args([
            "--mmseqs-db", "db"])
    end
end
