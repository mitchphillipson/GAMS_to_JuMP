$title	Grade Inflation -- L1 Estimation

set	seq /Overall/;

set	s	Students /
	Alex, Andy, Ariel, Billy, Bobby, Brett, Brook, Cameron, Cary, Casey,
	Chris, Dale, Dara, Darcy, Daryl, Devyn, Dominique, Drew, Emerson, Esme,
	Harley, Jade, Jordan, Kelly, Kim, Lindsey, Lou, Max, Meryl, Morgan,
	Peyton, Porntip, Reese, Robin, Sam, Skye, Sunny, Sydney, Tanner, Tracy /;

set	c	Courses /
	Apology, Astrology, Chronology, Cosmology, Demonology, Ecology,
	Epistemology, Etymology, Eulogy, Genealogy, Geology, Gynecology,
	Ideology, Immunology, Methodology, Morphology, Nephrology, Ontology,
	Pathology, Pharmacology, Philology, Phrenology, Psychology,
	Scientology, Seismology, Sociology, Statistics, Tautology, Technology,
	Terminology, Theology, Topology, Urology /;

set	g	Grades /A+, A, A-, B+, B, B-, C+, C, C-, D+, D, D-, F/;

set	marks(s,c,g) /
$ondelim
$include grades.txt
$offdelim
/;

parameter mark(g)	Translation of grades to marks /	
				A+	4.3
				A	4.0
				A-	3.7
				B+	3.3
				B	3.0
				B-	2.7
				C+	2.3
				C	2.0
				C-	1.7
				D+	1.3
				D	1
				D-	0.7
				F	0 /;

parameter	grade(s,c)	Numeric mark;

grade(s,c) = sum(marks(s,c,g), mark(g));


set	r(s,c)		Indication -- registration of student s in course c;

r(s,c) = yes$sum(marks(s,c,g),1);

alias (i,s), (j,c);

NONNEGATIVE
VARIABLE		T(i,j)	Absolute deviation;

VARIABLE		E(j)	Easiness of course j
			A(i)	Aptitude of student i,
			Z	Objective (least-squares calibration);


equations objdef, aveeasy, posdev, negdev;

objdef..		Z =g= sum((i,j),T(i,j));

aveeasy..		sum(j, E(j)) =e= 0;

posdev(r(i,j))..	T(i,j) =g= grade(i,j) - (A(i) + E(j));

negdev(r(i,j))..	T(i,j) =g= A(i) + E(j) - grade(i,j);

model grading /all/;
solve grading using LP minimizing Z;


display A.L, E.L;

parameter	results	Summary of results;
results(s,"Overall","actual") = sum(r(s,c), grade(s,c))/sum(r(s,c), 1);
results(s,"Overall","model") = A.L(s);
results(r(s,c),"actual") = grade(s,c);
results(r(s,c),"model") = A.L(s) + E.L(c);
option results:1;
display results;
