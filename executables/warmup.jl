function warmup_case(main, path)
    println("Warming up ", path)
    mktempdir() do output
        status = main([path, "--number-of-steps", "1", "--output-path", output, "--info-level", "-1"])
        iszero(status) || error("Warmup for $path exited with status $status")
    end
    return nothing
end
