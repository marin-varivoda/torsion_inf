AttachSpec("../third_party/mdmagma/v2/mdmagma.spec");
SetSeed(1337);

N := 53;
p := 2;
C := MDX1(N, GF(p));
FC := FunctionField(Curve(C));
g := 92; // genus of X1(53)

plc1  := Places(Curve(C),  1);
plc9  := Places(Curve(C),  9);

plc1sum := &+ plc1;

// [ 26, 0, 0, 0, 0, 0, 0, 0, 78]
print("The number of places of X1(53)(F2) of degrees 1..9 is:");
[#plc1] cat [#Places(Curve(C), i) : i in [2..8]] cat [#plc9];


plc9_diamond_orbits := [];
rep9 := [];
for deg9_place in plc9 do
    already_accounted_for := false;

    for diamond_orbit in plc9_diamond_orbits do
        if deg9_place in diamond_orbit then
            already_accounted_for := true;
            break;
        end if;
    end for;

    if not already_accounted_for then
        Append(~plc9_diamond_orbits, DiamondOrbit(C, deg9_place));
        Append(~rep9, deg9_place);
    end if;
end for;

// Degree 9 places partition into 3 diamond orbits. Their sizes: [ 26, 26, 26 ]
printf "Degree 9 places partition into %o diamond orbits. Their sizes: %o\n", #plc9_diamond_orbits, [#diamond_orbit : diamond_orbit in plc9_diamond_orbits];


// S := { D : D >=0, deg D = 18, and supp D has at least 9 different divisors of degree 1 }
// We build a list of divisors which dominates S up to diamond operator action.

// This means that for every s \in S, there should exist t \in divisors_to_check and a 
// diamond operator <a> such that <a> s <= t
divisors_to_check := [];


// 1) Cover divisors in S which have a place of degree 9 in support

// Divisor of this class has the form D = A + P where P is degree 9 place. Since D is required
// to have at least 9 different place 1 divisors in support, we observe that A is a sum of 
// exactly 9 degree 1 places, all with multiplicity 1. Therefore D <= plc1sum + P

for rep_place in rep9 do
    D := plc1sum + rep_place;
    Append(~divisors_to_check, D);
end for;


// 2) Cover divisors in S which are only supported on degree 1 places

// Since plc1 is entirely made up of cusps and diamond operators act transitively on this set,
// we can WLOG assume that D \in S has max multiplicity at plc1[1]. 

// Given D \in S with max multiplicity at plc1[1], we can write it as
// 
//      D = sum_{i=1}^{26} m_i*P_i
// 
// Define the "excess divisor" as 
// 
//      E = sum_{i=1}^{26} e_i*P_i, where e_i = max(m_i - 3, 0)
//
// If h = #supp E, we have: 18 = deg D >= 4*h + (9-h), which implies h <= 3

base_div := 3 * plc1sum; 

// 2.1) Excess vector E has 3 places in support, i.e., h = 3

// E must be of the form E = P_1 + P_i + P_j
// Therefore, we have D <= 3*plc1sum + plc1[1] + plc1[i] + plc1[j], for some i, j
for i in [2..26] do
    for j in [i+1 .. 26] do
        Append(~divisors_to_check, base_div + plc1[1] + plc1[i] + plc1[j]);
    end for;
end for;


// 2.2) h = 2

// E is of the form E = a*P_1 + b*P_i, where a >= b >= 1, and, (a+3)+(b+3)+7 <= 18, therefore a+b <= 5
// It follows that a<=4 and b<=2
// So, D <= 3*plc1sum + 4*plc1[1] + 2*plc1[i] for some i
for i in [2..26] do
    Append(~divisors_to_check, base_div + 4*plc1[1] + 2*plc1[i]);
end for;

// 2.3) h = 1

// E is of the form E = a * P_1 where (a+3)+8 <= 18, therefore a <= 7.
// So, D <= 3*plc1sum + 7*plc1[1]
Append(~divisors_to_check, base_div + 7*plc1[1]);

// 2.4) h=0

// This means that E=0, thus D <= 3*plc1sum, we don't need to add additional divisors to covering set




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

printf "Searching RR spaces for %o divisors... (this might take a moment)\n", #divisors_to_check;

T := Time();
results := ParallelMapRRSpaceHasFuncOfDegAtMost18(divisors_to_check);

num_true := #[res : res in results | res];
num_false := #results - num_true;

printf "Computation finished. %o divisors produced a non-constant function of degree <= 18 (true), and %o did not (false).\n", num_true, num_false;
printf "Total calculation time: %o seconds\n", Time(T);

exit;

