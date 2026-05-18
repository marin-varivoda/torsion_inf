// The following code computes the numerical values in the main theorem and Table 1 of the paper.


I:=Open("gonality1.dat", "r");
_,gonality1:=ReadObjectCheck(I);

I:=Open("gonality2.dat", "r");
_,gonality2:=ReadObjectCheck(I);

I:=Open("gonality3.dat", "r");
_,gonality3:=ReadObjectCheck(I);

I:=Open("bielliptic.dat", "r");
_,bielliptic:=ReadObjectCheck(I);

I:=Open("quintic.dat", "r");
_,quintic:=ReadObjectCheck(I);

assert &and [r`genus le 13: r in gonality1 cat gonality2 cat gonality3 cat bielliptic cat quintic];

assert #gonality1 eq 132;
assert #gonality2 eq 524;
assert #gonality3 eq 489; 
assert #bielliptic eq 1090;
assert #quintic eq 2;

assert [ #[r: r in gonality1 | r`genus eq g] : g in [0..14] ] eq [ 132, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
assert [ #[r: r in gonality2 | r`genus eq g] : g in [0..14] ] eq [ 0, 187, 177, 99, 12, 34, 2, 6, 1, 3, 0, 3, 0, 0, 0 ];
assert [ #[r: r in gonality3 | r`genus eq g] : g in [0..14] ] eq [ 0, 0, 0, 185, 249, 1, 24, 5, 16, 0, 8, 0, 1, 0, 0 ];
assert [ #[r: r in bielliptic | r`genus eq g] : g in [0..14] ] eq [ 0, 187, 132, 267, 173, 179, 21, 79, 5, 23, 18, 4, 0, 2, 0 ];
