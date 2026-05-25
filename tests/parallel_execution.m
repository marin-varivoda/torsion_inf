function isEven(x)
    return x mod 2 eq 0;
end function;

/**
 *  Parallel implementation for isEven. The idea is that following two calls
 *  produce the same result:
 *
 *      results := ParallelMapIsEven(task_inputs);
 *      results := [isEven(x) : x in task_inputs];
 *
 *  but ParallelMapIsEven performs the work in parallel.
 */
function ParallelMapIsEven(task_inputs)
    CORE_COUNT := 16;
    MEMORY_PER_WORKER := 16 * 10^9; // 0 for unlimited

    num_tasks := #task_inputs;
    if num_tasks eq 0 then
        return [];
    end if;

    worker_count := Minimum(CORE_COUNT, num_tasks);

    server_socket := Socket(: LocalHost := "localhost");
    host, port := Explode(SocketInformation(server_socket));

    // Start worker processes. Each worker receives a task index, computes the
    // corresponding result, and sends the pair "task_idx,result" back.
    for worker_idx in [1..worker_count] do
        child_pid := Fork();

        if child_pid eq 0 then
            SetMemoryLimit(MEMORY_PER_WORKER); 
            client_socket := Socket(host, port);
            
            // Ask for the first task. Index 0 means no completed task yet.
            Write(client_socket, "0,false");

            while true do
                _ := WaitForIO([client_socket]);
                msg := Read(client_socket);
                
                task_idx := -1;
                try
                    task_idx := eval msg;
                catch e
                    break;
                end try;

                // The parent sends -1 when there are no more tasks.
                if task_idx eq -1 then
                    break;
                end if;

                task_input := task_inputs[task_idx];
                task_result := isEven(task_input);
                result_str := Sprintf("%o", task_result);

                try
                    Write(client_socket, Sprintf("%o,%o", task_idx, result_str));
                catch e
                    break;
                end try;
            end while;

            // Terminate the child process to prevent normal cleanup on exit. Because fork()
            // means processes share open file descriptors, letting the child quit normally
            // would affect the file descriptors of the parent, possibly interfering with the interpreter.
            System(Sprintf("kill -9 %o", Getpid()));
            quit;
 
        end if;
    end for;

    worker_sockets := [];
    for worker_idx in [1..worker_count] do
        Append(~worker_sockets, WaitForConnection(server_socket));
    end for;

    results_by_task_idx := AssociativeArray();
    next_task_idx := 1;
    tasks_completed := 0;

    start_time := Realtime();
    report_interval := Maximum(1, Floor(num_tasks / 200)); // roughly every 0.5%
    printf "Starting parallel computation with %o tasks across %o workers...\n", num_tasks, worker_count;

    while tasks_completed lt num_tasks do
        if #worker_sockets eq 0 then
            error "Parallel computation failed: no worker sockets remain.";
        end if;

        ready_sockets := WaitForIO(worker_sockets);

        for sock in ready_sockets do
            try
                msg := Read(sock);
            catch e
                error "Parallel computation failed while reading from a worker.";
            end try;
            
            task_idx := -1;
            task_result := false;

            try
                comma_idx := Index(msg, ",");
                if comma_idx le 1 or comma_idx ge #msg then
                    error "Invalid worker message.";
                end if;

                task_idx := eval msg[1..comma_idx-1];
                result_str := msg[comma_idx+1..#msg];

                if task_idx lt 0 or task_idx gt num_tasks then
                    error "Invalid task index.";
                end if;

                if result_str eq "true" then
                    task_result := true;
                elif result_str eq "false" then
                    task_result := false;
                else
                    error "Invalid task result.";
                end if;
            catch e
                error "Parallel computation failed: worker returned an invalid message.";
            end try;

            if task_idx gt 0 then
                results_by_task_idx[task_idx] := task_result;
                tasks_completed +:= 1;

                if tasks_completed mod report_interval eq 0 or tasks_completed eq num_tasks then
                    elapsed := Realtime() - start_time;
                    core_time_per_task := (elapsed * worker_count) / tasks_completed;
                    eta := (num_tasks - tasks_completed) * (elapsed / tasks_completed);
                    pct := (tasks_completed * 100.0) / num_tasks;
                    
                    // We round the duration values to int so the output fits on the same line
                    printf "%o/%o (%o%%) | elapsed=%os | ETA=%os | core-time/task=%os\n",
                            tasks_completed, num_tasks, RealField(4)!pct,
                            Floor(elapsed), Floor(eta), Floor(core_time_per_task);
                end if;
            end if;

            if next_task_idx le num_tasks then
                try
                    Write(sock, Sprintf("%o", next_task_idx));
                catch e
                    error "Parallel computation failed while sending a task to a worker.";
                end try;
                next_task_idx +:= 1;
            else
                try
                    Write(sock, "-1");
                catch e
                    error "Parallel computation failed while stopping a worker.";
                end try;

                // This worker has no more tasks to receive. Remove its socket from the active set
                Exclude(~worker_sockets, sock);
            end if;
        end for;
    end while;

    WaitForAllChildren();
    
    return [ results_by_task_idx[i] : i in [1..num_tasks] ];
end function;

task_inputs := [1..1000000];

printf "Testing ParallelMapIsEven on %o inputs...\n", #task_inputs;

T := Time();
results := ParallelMapIsEven(task_inputs);
expected_results := [isEven(x) : x in task_inputs];

assert #results eq #task_inputs;
assert results eq expected_results;

num_true := #[res : res in results | res];
num_false := #results - num_true;

printf "Test passed. %o even numbers and %o odd numbers.\n", num_true, num_false;
printf "Total calculation time: %o seconds\n", Time(T);

exit;
