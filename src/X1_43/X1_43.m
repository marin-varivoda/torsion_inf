AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");
SetSeed(1337);

N := 43;
p := 2;
C := MDX1(N, GF(p));
FC := FunctionField(Curve(C));
g := 57; // genus of X1(43)

cusps := Cusps(C);
cusps_count := #cusps;
cusps_degree := [Degree(cusp) : cusp in cusps];

plc1  := Places(Curve(C),  1);
plc7  := Places(Curve(C),  7);
plc8  := Places(Curve(C),  8);
plc9  := Places(Curve(C),  9);
plc10 := Places(Curve(C), 10);
plc11 := Places(Curve(C), 11);

plc1sum := &+ plc1;

// Note that there are 21 cusps of degree 1 and 3 cusps of degree 7
// [ 21, 0, 0, 0, 0, 0, 6, 42, 105, 84, 126 ]
print("The number of places of X1(43)(F2) of degrees 1..11 is:");
[#plc1] cat [#Places(Curve(C), i) : i in [2..6]] cat [#plc7, #plc8, #plc9, #plc10, #plc11];


// S := { D : D >=0, deg D = 18, and supp D has at least 7 different divisors of degree 1 }
// We build a list of divisors which dominates S up to diamond operator action.

// This means that for every s \in S, there should exist t \in divisors_to_check and a 
// diamond operator <a> such that <a> s <= t
divisors_to_check := [];


// 1) Cover divisors in S which have a place of degree 10 or 11 in support 

/**
 *  Takes a list of places which is closed under diamond operator action,
 *  partitions it into orbits and returns a list where there is exactly one representative
 *  for each diamond orbit.
 */
function getDiamondOrbitRepresentatives(place_list)
    diamond_orbits := [];
    rep := [];
    for plc in place_list do
        already_accounted_for := false;

        // check if plc already is in some orbit
        for diamond_orbit in diamond_orbits do
            if plc in diamond_orbit then
                already_accounted_for := true;
                break;
            end if;
        end for;

        if not already_accounted_for then
            Append(~diamond_orbits, DiamondOrbit(C, plc));
            Append(~rep, plc);
        end if;
    end for;

    return rep;
end function;

/* 1.1) Cover divisors in S with degree 11 place in support */
rep11 := getDiamondOrbitRepresentatives(plc11); // 6 elements
for P in rep11 do
    D := plc1sum + P;
    Append(~divisors_to_check, D);
end for;

/* 1.2) Cover divisors in S with degree 10 place in support */
rep10 := getDiamondOrbitRepresentatives(plc10); // 4 elements
for P in rep10 do
    D := 2 * plc1sum + P;
    Append(~divisors_to_check, D);
end for;


// 2) Cover divisors in S which have a place of degree 8 or 9 in support

// Observe that a divisor of this class has the form A + P where P is degree 8 or 9 place
// and A >= 0 is supported on places of degree 1.

// By using diamond operators which act transitively on places of degree 1, we can
// WLOG assume that plc1[1] has highest multiplicity in A, counting argument shows that this
// multiplicity is at most 4. 

// Also it is impossible that A has two different places with multiplicity >= 3 since that would 
// mean the degree of the entire A + P is at least >= (3+3+1+1+1+1+1) + (8) = 19 which is impossible.  

// Therefore A + P <= 2*plc1sum + 2*plc1[1] + P.

// Finally, we store all the divisors which can appear on the right side of above inequality:

for P in plc9 do
    D := 2*plc1sum + 2*plc1[1] + P;
    Append(~divisors_to_check, D);
end for;

for P in plc8 do
    D := 2*plc1sum + 2*plc1[1] + P;
    Append(~divisors_to_check, D);
end for;


// 3) Cover divisors in S which have a place of degree 7 in support

// The form is A + P again where P is a degree 7 place and A >= 0 is supported on places of degree 1.
// WLOG assume that plc1[1] has highest multiplicity in A, counting argument shows that this
// multiplicity is at most 5. 

// Moreover, at most two places in A can have multiplicity >= 3, and in that case both have 
// exactly the multiplicity 3. Note that (3+3+1+1+1+1+1) + (7) = 18

// Therefore A + P <= 2*plc1sum + 3*plc1[1] + plc1[i] + plc7[j], for some i,j 

for i in [2..#plc1] do
    for j in [1..#plc7] do
        D := 2*plc1sum + 3*plc1[1] + plc1[i] + plc7[j];
        Append(~divisors_to_check, D);
    end for;
end for;


// 4) Cover divisors in S which only have degree 1 places in support

// We use a separate C++ file to generate the cover in this case.
// The file X1_43/case4_cover.txt has one divisor per line:
//      a1 a2 ... a21
// representing a1*plc1[1] + ... + a21*plc1[21].

function LoadCase4Cover(filename, plc1)
    assert #plc1 eq 21;

    FP := Open(filename, "r");
    ret := [];

    while true do
        line := Gets(FP);

        if IsEof(line) then
            break;
        end if;

        // Skip empty lines, just in case.
        if #line eq 0 then
            continue;
        end if;

        coeffs := StringToIntegerSequence(line);

        assert #coeffs eq #plc1;
        assert &and[c ge 0 : c in coeffs];

        D := &+[ coeffs[i] * plc1[i] : i in [1..#plc1] ];
        assert Degree(D) le 57;

        Append(~ret, D);
    end while;

    return ret;
end function;

case4_cover := LoadCase4Cover("X1_43/case4_cover.txt", plc1);
assert #case4_cover eq 694; // make sure the parsed size matches C++ output

printf "Loaded case 4 cover of size: %o\n", #case4_cover; 

divisors_to_check := divisors_to_check cat case4_cover;

/**
 *  Returns true/ false depending on whether or not RR Space has a non-constant function
 *  of degree at most 18
 */
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

/**
 *  Parallel implementation for RRSpaceHasFuncOfDegAtMost18. The idea is that 
 *  following two calls produce the same result:
 *  
 *  results := ParallelMapRRSpaceHasFuncOfDegAtMost18(divisors_to_check);
 *  results := [RRSpaceHasFuncOfDegAtMost18(D) : D in divisors_to_check];
 *  
 *  but ParallelMapRRSpaceHasFuncOfDegAtMost18 will do this work in parallel
 *  therefore being faster in practice on multicore systems.
 */
function ParallelMapRRSpaceHasFuncOfDegAtMost18(task_inputs)
    CORE_COUNT := 32;
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
                task_result := RRSpaceHasFuncOfDegAtMost18(task_input);
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
                    
                    printf "Progress: %o%% (%o/%o) | Elapsed: %os | ETA: %os | Avg/Task/Core: %os\n", 
                        RealField(4)!pct, tasks_completed, num_tasks, 
                        RealField(5)!elapsed, RealField(5)!eta, RealField(5)!core_time_per_task;
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

printf "Searching RR spaces for %o divisors... (this might take a moment)\n", #divisors_to_check;

T := Time();
results := ParallelMapRRSpaceHasFuncOfDegAtMost18(divisors_to_check);

num_true := #[res : res in results | res];
num_false := #results - num_true;

printf "Computation finished. %o divisors produced a non-constant function of degree <= 18 (true), and %o did not (false).\n", num_true, num_false;
printf "Total calculation time: %o seconds\n", Time(T);

exit;