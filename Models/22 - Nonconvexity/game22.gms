$title	Board Game "22" -- Slow Version Ignores Symmetry

set	n0	All nodes (include center node 0) /n0*n12/,
	n(n0)	Nodes other than node 0 /n1*n12/
	p	Pieces /1*13/,
	t	Triangles /t1*t6/

	tri(t,n) Assignment of triangles to nodes
		 / t1.(n1,n3,n4),  t2.(n4,n5,n7), t3.(n7,n10,n11),
		   t4.(n9,n10,n12), t5.(n6,n8,n9), t6.(n2,n3,n6)/,

	s	All possible solutions /sol1*sol5000/,

	ss(s)	Solutions which have been found, 

	cut(s,p,n)	Cut which makes previously detected solutions infeasible;

alias (p,pp);

variable	OBJ	Vacuous objective;

binary
variable	Z(p,n)	Assignment of piece p to node n;

equations	objdef, used, filled, linesum, othersol;

objdef..	OBJ =e= 0;

used(p)..	sum(n, Z(p,n)) =L= 1;

filled(n)..	sum(p, Z(p,n)) =E= 1;

linesum(t)..	sum((p,tri(t,n)), p.val*Z(p,n)) =e= 22;

othersol(ss)..	sum(cut(ss,p,n), Z(p,n)) =l= 11;

cut(s,p,n) = no;
ss(s) = no;

model ip /all/;
solve ip using mip maximizing obj;

parameter	nsol		Solution count /0/,
		sollist(p,*,*)	Solution assignments;

*	Restrict model output to improve performance:

ip.limrow=0; ip.limcol=0; option solprint=off; 

*	We need to solve a number of models -- to make this go
*	more rapidly, let GAMS load the solver in memory to 
*	avoid the cost of reloading the library each time:

ip.solvelink = %solvelink.LoadLibrary%;

*	Outerloop over sol1(s) provides an "anchor point" for
*	assignments to the offet reference s+nsol:

loop(s$sameas(s,"sol1"),

*	Find solutions with peg pp in the center hole:

  loop(pp,

*	Omit pp from the assignment into nodes other than 0:

	Z.FX(pp,n) = 0;

*	First solution:

	solve ip using mip maximizing obj;

*	While we have a solution and we have found 9 or 
*	fewer solution with peg pp in the center, continue
*	solving:

	while (((ip.modelstat<>10) and (nsol<=9)),


*	We have a solution in hand at this.  Add this to the 
*	list of cuts, using "s+nsol" to refer to this new solution:

	  cut(s+nsol,p,n)$Z.L(p,n) = yes;

*	Add the new solution to the list of cuts:

	  ss(s+nsol) = yes;

*	Keep record of which pegs are in which holes for this
*	solution:

	  sollist(pp,"n0",s+nsol) = pp.val;
	  loop((p,n)$Z.L(p,n), sollist(pp,n,s+nsol) = p.val; );

*	See if we can find another solution:

	  nsol = nsol + 1;
	  solve ip using mip maximizing obj;

	);

*	Done with peg pp in hole 0.  Reset the solution count
*	when we move to putting the next peg in hole 0:

	Z.UP(pp,n) = 1;
	ss(s) = no;
	nsol = 0;
));

option sollist:0:2:1;
display sollist;

parameter	solchk	Cross check on the solution;
solchk(p,s,t)$sum(tri(t,n), sollist(p,n,s)) = sum(tri(t,n), sollist(p,n,s)) - 22;
option solchk:0:0:1;
display solchk;

$exit

execute_unload 'game22.gdx',sollist;
execute 'gdxxrw i=game22.gdx o=game22.xlsx par=sollist rng=Solutions rdim=1 cdim=1 clear';
