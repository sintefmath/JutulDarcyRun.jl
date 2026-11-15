function tee_run_stream(pipe, terminal, log)
    try
        while !eof(pipe.out)
            bytes = readavailable(pipe.out)
            write(log, bytes)
            flush(log)
            try
                write(terminal, bytes)
                flush(terminal)
            catch err
                # Continue writing the log if a terminal pipe closes early.
                err isa Base.IOError || rethrow()
                terminal = devnull
            end
        end
    finally
        close(pipe.out)
    end
    return nothing
end

function with_run_log(run, output_dir, name, input)
    mkpath(output_dir)
    path = joinpath(output_dir, name * ".DBG")
    return mkpidlock(path * ".lock"; stale_age = 60) do
        open(path, "w") do log
            println(log, "Running ", input)
            flush(log)
            out, err = Pipe(), Pipe()
            Base.link_pipe!(out; reader_supports_async = true, writer_supports_async = true)
            Base.link_pipe!(err; reader_supports_async = true, writer_supports_async = true)
            out_terminal, err_terminal = stdout, stderr
            out_task = @async tee_run_stream(out, out_terminal, log)
            err_task = @async tee_run_stream(err, err_terminal, log)
            try
                result = try
                    redirect_stdio(; stdout = out, stderr = err) do
                        try
                            run()
                        finally
                            Base.Libc.flush_cstdio()
                        end
                    end
                finally
                    closewrite(out)
                    closewrite(err)
                    wait(out_task)
                    wait(err_task)
                end
                println(log, "Run complete")
                return result
            catch exception
                print(log, "Run failed. ")
                showerror(log, exception)
                println(log)
                rethrow()
            end
        end
    end
end
