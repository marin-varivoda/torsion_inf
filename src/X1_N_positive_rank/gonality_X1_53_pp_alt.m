A<x,y> := PolynomialRing(GF(2),2);
P,homogenize := Homogenization(A);
print("f := The polynomial from math.mit.edu/~drew/FFFc53.txt reduced mod 2");
f := y^12*(y+1)^4*x^39+y^12*(y+1)^5*x^38+y^6*(y^18+y^15+y^11+y^10+y^5+y^4+y^3+y^2+1)*(y^4+y^3+y^2+y+1)*(y+1)^3*x^25+y^5*(y^10+y^9+y^8+y^6+y^2+y+1)*(y^15+y^14+y^12+y^10+y^9+y^5+y^2+y+1)*(y+1)^3*x^23+(y^3+y^2+1)*(y^23+y^22+y^20+y^18+y^16+y^15+y^9+y^8+y^5+y^4+1)*y^5*(y+1)^3*x^22+(y^3+y^2+1)^2*y^11*(y^4+y^3+1)*(y^2+y+1)*x^33+(y^3+y^2+1)^2*(y^8+y^7+y^5+y+1)*y^8*(y+1)^3*x^31+y^7*(y^5+y^4+y^3+y^2+1)*(y^11+y^9+y^7+y^4+1)*(y+1)^5*x^28+y^8*(y^5+y^3+y^2+y+1)*(y^7+y^6+y^5+y^4+y^2+y+1)*(y+1)^9*x^27+y^11*(y^6+y^5+y^4+y+1)*x^36+y^5*(y^31+y^29+y^28+y^25+y^22+y^21+y^20+y^19+y^18+y^17+y^16+y^13+y^12+y^6+y^2+y+1)*x^20+(y^11+y^9+y^7+y^6+y^4+y^2+1)*(y+1)^14*x^5+y^5*(y+1)^16*x+(y^5+y^3+1)*y^8*(y^2+y+1)*(y^7+y^6+y^4+y^2+1)*(y+1)^8*x^26+y*(y^5+y^4+y^3+y^2+1)^2*(y^13+y^7+y^5+y^3+y^2+y+1)*(y+1)*(y^9+y^8+y^7+y^3+y^2+y+1)*x^16+y*(y^22+y^21+y^20+y^19+y^18+y^17+y^11+y^10+y^9+y^8+y^5+y^2+1)*(y^2+y+1)*(y^6+y^3+1)*(y+1)^4*x^15+(y^3+y+1)*y*(y^17+y^16+y^14+y^12+y^7+y^5+y^4+y+1)*(y^7+y^6+y^5+y^4+1)*(y+1)^5*x^14+y*(y^15+y^14+y^11+y^6+1)*(y+1)^6*(y^5+y^4+y^3+y+1)*(y^6+y^5+y^4+y^2+1)*x^13+y*(y^5+y^2+1)*(y^7+y^5+y^4+y^3+y^2+y+1)*(y+1)^9*(y^14+y^9+y^8+y^7+y^5+y^3+y^2+y+1)*x^10+y^5*(y^2+y+1)*(y+1)^22*x^3+(y^15+y^14+y^13+y^12+y^10+y^8+y^7+y^5+y^4+y^2+1)*y^8*(y+1)^3*x^30+(y^8+y^6+y^5+y^4+y^2+y+1)*y^7*(y+1)^12*x^29+y^5*(y^24+y^23+y^22+y^21+y^19+y^18+y^16+y^15+y^12+y^11+y^10+y^9+y^8+y^6+y^5+y^4+1)*(y+1)^3*x^24+(y^18+y^16+y^15+y^14+y^12+y^11+y^8+y^5+y^2+y+1)*y^4*(y^14+y^13+y^12+y^11+y^10+y^9+y^8+y^6+y^5+y^4+1)*x^18+(y^25+y^24+y^23+y^22+y^19+y^18+y^17+y^16+y^15+y^14+y^12+y^11+y^4+y+1)*y^2*(y^8+y^7+y^6+y^4+y^2+y+1)*x^17+y^11*(y^2+y+1)*(y+1)^4*x^37+(y^3+y^2+1)*y^12*(y^5+y^4+y^3+y^2+1)*x^35+y^12*(y^5+y^4+y^3+y^2+1)*(y^5+y^4+y^2+y+1)*x^34+y^11*(y^12+y^11+y^8+y^2+1)*(y+1)*x^32+y^4*(y+1)^17+y^5*(y^5+y^4+y^3+y^2+1)*(y^21+y^20+y^17+y^15+y^14+y^12+y^9+y^8+y^7+y^6+y^4+y^2+1)*(y+1)^4*x^21+y^4*(y^15+y^12+y^11+y^8+y^6+y^5+y^4+y^3+y^2+y+1)*(y^6+y^3+1)*(y^12+y^11+y^10+y^9+y^7+y^5+y^4+y^3+1)*x^19+(y^12+y^11+y^10+y^7+y^6+y^4+y^2+y+1)*y*(y^16+y^14+y^9+y^6+y^3+y+1)*(y+1)^7*x^12+(y^12+y^11+y^8+y^7+y^6+y^5+y^4+y^3+1)*y*(y^16+y^13+y^12+y^11+y^10+y^8+y^5+y+1)*(y+1)^8*x^11+(y^5+y^3+1)*(y^13+y^12+y^11+y^9+y^7+y^5+y^4+y^3+1)*y*(y+1)^10*x^9+(y^3+y^2+1)*(y^14+y^12+y^10+y^8+y^6+y^5+y^4+y^2+1)*y*(y+1)^11*x^8+(y^14+y^13+y^11+y^8+y^6+y^5+y^4+y^2+1)*y*(y^2+y+1)*(y+1)^12*x^7+y^2*(y^6+y^4+y^3+y+1)*(y^4+y+1)*(y+1)^13*x^6+(y^3+y^2+1)*y^4*(y^6+y+1)*(y+1)^15*x^4+(y^3+y^2+1)*(y^3+y+1)*y^6*(y+1)^16*x^2 ;

f := homogenize(f);
XF2 := Scheme(ProjectiveSpace(P),f);
FF := FunctionField(XF2);
AFF := AlgorithmicFunctionField(FF);
plc1 := Places(AFF,1);
cuspsum := &+ plc1;
plc9 := Places(AFF, 9);

divisors_to_check := [];

print "Adding divisors consisting of degree 1 places.";


base_div := 3 * cuspsum; 

// Covers the case of a single heavy cusp (max pole up to 10)
Append(~divisors_to_check, base_div + 7*plc1[1]);
print "Added Profile A. Total so far:", #divisors_to_check;


// Profile B: E1 = 4, E2 = 2 (The Super-Profile for 2 heavy cusps) -> 25 combinations
for i in [2..26] do
    Append(~divisors_to_check, base_div + 4*plc1[1] + 2*plc1[i]);
end for;
print "Added Profile B. Total so far:", #divisors_to_check;

// Stage 3: Profile C (E1 = 1, E2 = 1, E3 = 1)
// Profile C: E1 = 1, E2 = 1, E3 = 1 (For 3 heavy cusps) -> 300 combinations
for i in [2..26] do
    for j in [i+1 .. 26] do
        Append(~divisors_to_check, base_div + 1*plc1[1] + 1*plc1[i] + 1*plc1[j]);
    end for;
end for;
print "Added Profile C. Final list size:", #divisors_to_check;
print "-----------------------------------------";

function RRSpaceHasFuncOfDegAtMost18(D)
    R, m := RiemannRochSpace(D);
    for f in R do
        g := m(f);
        deg_f := Degree(g);
        if deg_f gt 0 and deg_f le 18 then
            return true;
        end if;
    end for;
    return false;
end function;

print "Doing a mock Riemann-Roch computation to build cache before running workers.";
tmp := RRSpaceHasFuncOfDegAtMost18(cuspsum);
print "Mock computation finished.";

print "Initiating parallel computation!!!";

// ======================================================================
// BEGIN AUTO-GENERATED PARALLEL BOILERPLATE (MACRO ID: 1)
// ======================================================================
function _RunParallelSocketMap_1(task_inputs_list)
    CORE_COUNT := 24; // Auto-detected by Python preprocessor
    num_tasks := #task_inputs_list;

    server_socket := Socket(: LocalHost := "localhost");
    host, port := Explode(SocketInformation(server_socket));

    // --- SPAWN WORKERS ---
    for i in [1..CORE_COUNT] do
        child_pid := Fork();

        if child_pid eq 0 then
            SetMemoryLimit(16 * 10^9); 
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
                task_result := RRSpaceHasFuncOfDegAtMost18(input_val);
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

results := _RunParallelSocketMap_1(divisors_to_check);

print "Results: ", results;
