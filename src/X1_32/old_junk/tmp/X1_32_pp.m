AttachSpec("/home/mvarivoda/mdmagma/v2/mdmagma.spec");
SetSeed(1337);


N := 32;

C := MDX1(N, Rationals());
C_proj := Curve(C);

// Reduce the curve over finite field
p := 3;
C_fin := Curve(Reduction(C_proj, p));

printf "Initialized the model for X_1(%o) and reduced it over a finite field with %o elements.\n", N, p;
printf "Checking if reduced X_1(%o) has a function of degree at most 9... (this might take a moment)\n", N;

plcs := [];
for i in [1..9] do
    printf "Starting Places(Curve(C), %o) [%o/12]...\n", i, i;
    t0 := Cputime();

    Append(~plcs, Places(Curve(C_fin), i));

    printf "Finished Places(Curve(C), %o) [%o/12] in %o seconds. Found %o places.\n",
        i, i, Cputime(t0), #plcs[#plcs];
end for;


Div := DivisorGroup(Curve(C_fin));
T_divisors := [];

filename := "covering_set_T.txt";
printf "\nOpening %o to load the dominating set T...\n", filename;
t0 := Cputime();

// Open file for reading
F := Open(filename, "r");

count := 0;
line := Gets(F);

while not IsEof(line) do
    if #line gt 0 then
        tokens := Split(line, ",");
        
        // Build the divisor by summing the places
        d := Div ! 0;
        for token in tokens do
            // Parse strings like "P_9_1976" -> parts = ["P", "9", "1976"]
            parts := Split(token, "_");
            deg := StringToInteger(parts[2]);
            idx := StringToInteger(parts[3]);
            
            d +:= Div ! plcs[deg][idx];
        end for;
        
        Append(~T_divisors, d);
    end if;
    
    count +:= 1;
    if count mod 5000 eq 0 then
        printf "  Parsed %o divisors...\n", count;
    end if;
    
    line := Gets(F);
end while;

delete F; // Close the file

printf "Successfully loaded %o divisors into Magma in %o seconds.\n", #T_divisors, Cputime(t0);


function RRSpaceHasFuncOfDegExactly9(D)
    R, m := RiemannRochSpace(D);
    for f in R do
        g := m(f);
        deg_f := Degree(g);
        if deg_f eq 9 then
            return true;
        end if;
    end for;
    return false;
end function;

// ======================================================================
// BEGIN AUTO-GENERATED PARALLEL BOILERPLATE (MACRO ID: 1)
// ======================================================================
function _RunParallelSocketMap_1(task_inputs_list)
    CORE_COUNT := 30; // Auto-detected by Python preprocessor
    num_tasks := #task_inputs_list;

    server_socket := Socket(: LocalHost := "localhost");
    host, port := Explode(SocketInformation(server_socket));

    // --- SPAWN WORKERS ---
    for i in [1..CORE_COUNT] do
        child_pid := Fork();

        if child_pid eq 0 then
            SetMemoryLimit(2 * 10^9); 
            client_socket := Socket(host, port);
            
            // Ping the master to ask for the first task
            Write(client_socket, "0,false");

            while true do
                _ := WaitForIO([client_socket]);
                msg := Read(client_socket);
                
                // BULLETPROOFING: If master dies, eval will throw. Catch and exit.
                task_idx := -1;
                try
                    task_idx := eval msg;
                catch e
                    break;
                end try;

                // -1 is the signal from the master to shut down
                if task_idx eq -1 then
                    break;
                end if;

                // Do the work
                input_val := task_inputs_list[task_idx];
                task_result := RRSpaceHasFuncOfDegExactly9(input_val);
                serialized_result := Sprintf("%o", task_result);

                // Report back: "task_index,result"
                try
                    Write(client_socket, Sprintf("%o,%o", task_idx, serialized_result));
                catch e
                    break;
                end try;
            end while;

            // Failsafe teardown
            System(Sprintf("kill -9 %o", Getpid()));
            quit; 
        end if;
    end for;

    // --- MASTER EVENT LOOP ---
    worker_sockets := [];
    for i in [1..CORE_COUNT] do
        Append(~worker_sockets, WaitForConnection(server_socket));
    end for;

    results_assoc := AssociativeArray();
    next_task_to_assign := 1;
    tasks_completed := 0;

    // --- PROGRESS MONITORING SETUP ---
    t0 := Realtime();
    report_interval := Maximum(1, Floor(num_tasks / 200)); // 0.5% threshold
    printf "Starting macro %o with %o tasks across %o workers...\n", 1, num_tasks, CORE_COUNT;

    while tasks_completed lt num_tasks do
        ready_sockets := WaitForIO(worker_sockets);

        for sock in ready_sockets do
            msg := Read(sock);
            
            // --- THE BULLETPROOF FIX ---
            idx := -1;
            res := false;
            parse_success := false;
            try
                // Use Index instead of Split in case the serialized string contains commas
                comma_idx := Index(msg, ",");
                if comma_idx gt 0 then
                    idx := eval msg[1..comma_idx-1];
                    res_str := msg[comma_idx+1..#msg];
                    res := eval res_str;
                    parse_success := true;
                end if;
            catch e
                Exclude(~worker_sockets, sock); // Socket is dead
                continue;
            end try;
            // ---------------------------

            if parse_success then
                if idx gt 0 then
                    results_assoc[idx] := res;
                    tasks_completed +:= 1;

                    // --- PROGRESS REPORTING ---
                    if tasks_completed mod report_interval eq 0 or tasks_completed eq num_tasks then
                        elapsed := Realtime() - t0;
                        
                        // Time it takes for a SINGLE core to do one task
                        core_time_per_task := (elapsed * CORE_COUNT) / tasks_completed;
                        
                        // ETA based on overall system throughput
                        eta := (num_tasks - tasks_completed) * (elapsed / tasks_completed);
                        
                        pct := (tasks_completed * 100.0) / num_tasks;
                        
                        printf "Progress: %o%% (%o/%o) | Elapsed: %os | ETA: %os | Avg/Task/Core: %os\n", 
                            RealField(4)!pct, tasks_completed, num_tasks, 
                            RealField(5)!elapsed, RealField(5)!eta, RealField(5)!core_time_per_task;
                    end if;
                end if;
            end if;

            // Give the worker a new task OR tell it to quit
            if next_task_to_assign le num_tasks then
                try
                    Write(sock, Sprintf("%o", next_task_to_assign));
                catch e
                    Exclude(~worker_sockets, sock);
                end try;
                next_task_to_assign +:= 1;
            else
                try
                    Write(sock, "-1"); // No more tasks, signal worker to retire
                catch e
                    Exclude(~worker_sockets, sock);
                end try;
            end if;
        end for;
    end while;

    WaitForAllChildren();
    
    // Flatten the associative array into a standard Magma sequence
    return [ results_assoc[i] : i in [1..num_tasks] ];
end function;
// ======================================================================
// END AUTO-GENERATED PARALLEL BOILERPLATE (MACRO ID: 1)
// ======================================================================

results := _RunParallelSocketMap_1(T_divisors);

print results;

for r in results do
    if r then
        print "Found function of degree 9 :/";
    end if;
end for;
