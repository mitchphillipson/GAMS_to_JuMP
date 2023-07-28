$title	Who owns the fish?

$ontext

This problem was passed to my by my colleague Phil Graves.  I also need to 
thank Erwin Kalvalagen for teaching me what is meant by a "cut" in integer
programming.  

Thomas Rutherford
Economics Department
University of Colorado

GIVEN

1. In a street there are five houses, painted five different colours. 
2. In each house lives a person of different nationality 
3. These five homeowners each drink a different kind of beverage, smoke 
different brand of cigar and keep a different pet. 

THE QUESTION: WHO OWNS THE FISH? 

HINTS 

1. The Brit lives in a red house. 
2. The Swede keeps dogs as pets. 
3. The Dane drinks tea. 
4. The Green house is just to the left of the White house. 
5. The owner of the Green house drinks coffee. 
6. The person who smokes Pall Mall rears birds. 
7. The owner of the Yellow house smokes Dunhill. 
8. The man living in the centre house drinks milk. 
9. The Norwegian lives in the first house. 
10. The man who smokes Blends lives next to the one who keeps cats. 
11. The man who keeps horses lives next to the man who smokes Dunhill. 
12. The man who smokes Blue Master drinks beer. 
13. The German smokes Prince. 
14. The Norwegian lives next to the blue house. 
15. The man who smokes Blends has a neighbour who drinks water. 

It has been asserted that Albert Einstein wrote this riddle early
during the 19th century. He said that 98% of the world population
would not be able to solve it. 

Erwin Kalvalagen pointed me to a web site which comments on Einstein
being the author of this riddle (supposedly early on in his career):

"Finally, you may want to notice that the young Albert Einstein
(1879-1955) could not possibly have authored the puzzle in this form:
The Pall Mall brand of cigarettes was introduced by Butler & Butler in
1899 (sold to American Tobacco in 1907 and Brown & Williamson in 1994)
and Alfred Dunhill was established in 1893 (starting to manufacture
pipes in 1907) when Einstein was still a young man. However, the Blue
Master brand was introduced by J. L.  Tiedemann in 1937, when Einstein
was 58!"

$offtext

*	Define five sets defining the five characteristics of
*	each individual.  The way the program is designed, the symbols
*	used to define characteristics must be distinct (e.g., you
*	cannot define an element of the smokes set S named "Red").

set	h	House /h1*h5/
	c	House colors /Red, Green, Yellow, Blue, White/,
	s	Smokes / Pall-Mall, Dunhill, Blends, Prince, Blue-Master/,
	b	beverages / Coffee, Milk, Beer, Water, Tea/,
	p	Pet	/Dogs, Birds, Cats, Horses, Fish/,
	n	Nationality /Brit, Swede, Dane, Norwegian, German/;


*	We need a second reference to the houses in order to express hint 4:

alias (h,hh);


variable	obj	Objective function (vacuous);

nonnegative
variable
		Z(h,c,s,b,p,n)	Choice variable is 1 if we have h-c-s-b-p-n living on street;

equations	housing, colors, smokes, beverages, pets, nations,
			hint4,hint10,hint11,hint14,hint15,objdef;

*	One persone in each house:

housing(h)..	sum((c,s,b,p,n), Z(h,c,s,b,p,n)) =e= 1;

*	Each of five colors appears:

colors(c)..	sum((h,s,b,p,n), Z(h,c,s,b,p,n)) =e= 1;

*	Each of five types of smokes:

smokes(s)..	sum((h,c,b,p,n), Z(h,c,s,b,p,n)) =e= 1;

*	Each of five beverages:

beverages(b)..	sum((h,c,s,p,n), Z(h,c,s,b,p,n)) =e= 1;

*	Each of five varieties of pets:

pets(p)..	sum((h,c,s,b,n), Z(h,c,s,b,p,n)) =e= 1;

*	Each of five nations:

nations(n)..	sum((h,c,s,b,p), Z(h,c,s,b,p,n)) =e= 1;

*	4. The Green house is on the left of the White house. 

hint4(h)..		sum((s,b,p,n), Z(h,"Green",s,b,p,n)) =L= 
		sum((s,b,p,n), Z(h+1,"White",s,b,p,n));

*	10. The man who smokes Blends lives next to the one who keeps Cats. 

hint10(h)..	sum((c,b,p,n),Z(h,c,"Blends",b,p,n))
			=L= sum((c,s,b,n),  Z(h+1,c,s,b,"Cats",n))
			  + sum((c,s,b,n),  Z(h-1,c,s,b,"Cats",n));

*	11. The mann who keeps horses lives next to the man who smokes Dunhill. 

hint11(h)..	sum((c,s,b,n),Z(h,c,s,b,"Horses",n))
		 =L= sum((c,b,p,n),  Z(h+1,c,"Dunhill",b,p,n))
		   + sum((c,b,p,n),  Z(h-1,c,"Dunhill",b,p,n));

*	14. The Norwegian lives next to the blue house. 

hint14(h)..	sum((c,s,b,p),Z(h,c,s,b,p,"Norwegian"))
		 =L= sum((s,b,p,n),  Z(h-1,"Blue",s,b,p,n))
		   + sum((s,b,p,n),  Z(h+1,"Blue",s,b,p,n));

*	15. The man who smokes Blends has a neighbour who drinks Water. 

hint15(h)..	sum((c,b,p,n),Z(h,c,"Blends",b,p,n))
		 =L= sum((c,s,p,n),  Z(h-1,c,s,"Water",p,n))
		   + sum((c,s,p,n),  Z(h+1,c,s,"Water",p,n));

*	The objective function is not meaningful but required for the mip problem type:

objdef..	obj =e= 0;

model einstein /all/;

*	Then start to eliminate options based on those hints which link specific
*	pairs of traits.  If a hint indicates that two traits are linked, we can drop all
*	permutations in which one or the other but not both characteristics are present.

*	Z=1 makes an assignment:

Z.UP(h,c,s,b,p,n) = 1;

$ontext

Paul van der Eijk from GAMS writes:

"I was looking at the construct you used, like:

set	match(*,*)		Used to set bounds on Z;

match(c,n) = yes; 
match("Red",n) = no; 
match(c,"Brit")=no;  
match("Red","Brit") = yes;
z.fx(h,c,s,b,p,n)$(not match(c,n)) = 0;
match(c,n) = no;
...

I think this can be replaced by:

z.fx(h,c,s,b,p,n)$(SameAs(n,"Brit") xor SameAs(c,"Red")) = 0;

First time I see a good use of the xor operator!"

Bruce McCarl's GAMS user's guide provides the following description of the
XOR operator:

When one wishes to perform an action if and only if one of two or more
conditionals apply one can join them with an xor operator. This
involves using syntax like

Action$(logical condition 1 xor logical condition 2)

or

If(logical condition 1 xor logical condition 2 , Action);

or

While(logical condition 1 xor logical condition 2, Action);

$offtext

*	1. The Brit lives in a red house. 

Z.fx(h,c,s,b,p,n)$(SameAs(n,"Brit") xor SameAs(c,"Red")) = 0;

*	2. The Swede keeps Dogs as pets. 

Z.fx(h,c,s,b,p,n)$(SameAs(n,"Swede") xor SameAs(p,"Dogs")) = 0;

*	3. The Dane drinks tea. 

Z.fx(h,c,s,b,p,n)$(SameAs(n,"Dane") xor SameAs(b,"Tea")) = 0;


*	5. The owner of the Green house drinks coffee. 

Z.fx(h,c,s,b,p,n)$(SameAs(c,"Green") xor SameAs(b,"Coffee")) = 0;

*	6. The person who smokes Pall Mall rears Birds. 

Z.fx(h,c,s,b,p,n)$(SameAs(s,"Pall-Mall") xor SameAs(p,"Birds")) = 0;

*	7. The owner of the Yellow house smokes Dunhill. 

Z.fx(h,c,s,b,p,n)$(SameAs(c,"Yellow") xor SameAs(s,"Dunhill")) = 0;

*	8. The man living in the centre house drinks milk. 

Z.fx(h,c,s,b,p,n)$(SameAs(b,"Milk") xor SameAs(h,"H3")) = 0;

*	9. The Norwegian lives in the first house. 

Z.fx(h,c,s,b,p,n)$(SameAs(n,"Norwegian") xor SameAs(h,"H1")) = 0;

*	12. The man who smokes Blue Master drinks Beer. 

Z.fx(h,c,s,b,p,n)$(SameAs(b,"Beer") xor SameAs(s,"Blue-Master")) = 0;

*	13. The German smokes Prince. 

Z.fx(h,c,s,b,p,n)$(SameAs(n,"German") xor SameAs(s,"Prince")) = 0;

*	End of screen specific pairs.  We can now report the number of combinations
*	which remain:

scalar	npick		Number of combinations remaining from which to choose 5;
npick = sum((h,c,s,b,p,n)$Z.up(h,c,s,b,p,n), 1);
display npick;

*	Tell GAMS to omit variables from the model which are fixed.

einstein.holdfixed = yes;

solve einstein using lp minimizing obj;

option Z:0:0:1;
display Z.L;

$ontext

Exercise:

Revise Hint 4 to be:

4. The Green house is to the left of the White house. 

Modify the code to find all 7 integer solutions:

sol1	h1.Green .Pall-Mall  .Coffee.Birds .Norwegian
	h2.Blue  .Prince     .Water .Fish  .German   
	h3.Red   .Blends     .Milk  .Horses.Brit     
	h4.Yellow.Dunhill    .Tea   .Cats  .Dane     
	h5.White .Blue-Master.Beer  .Dogs  .Swede    
sol2	h1.Green .Pall-Mall  .Coffee.Birds .Norwegian
	h2.Blue  .Prince     .Water .Fish  .German   
	h3.White .Blends     .Milk  .Dogs  .Swede    
	h4.Yellow.Dunhill    .Tea   .Cats  .Dane     
	h5.Red   .Blue-Master.Beer  .Horses.Brit     
sol3	h1.Yellow.Dunhill    .Water .Cats  .Norwegian
	h2.Blue  .Blends     .Tea   .Horses.Dane     
	h3.Red   .Pall-Mall  .Milk  .Birds .Brit     
	h4.Green .Prince     .Coffee.Fish  .German   
	h5.White .Blue-Master.Beer  .Dogs  .Swede    
sol4	h1.Green .Pall-Mall  .Coffee.Birds .Norwegian
	h2.Blue  .Prince     .Water .Cats  .German   
	h3.Red   .Blends     .Milk  .Horses.Brit     
	h4.Yellow.Dunhill    .Tea   .Fish  .Dane     
	h5.White .Blue-Master.Beer  .Dogs  .Swede    
sol5	h1.Green .Pall-Mall  .Coffee.Birds .Norwegian
	h2.Blue  .Prince     .Water .Cats  .German   
	h3.White .Blends     .Milk  .Dogs  .Swede    
	h4.Yellow.Dunhill    .Tea   .Fish  .Dane     
	h5.Red   .Blue-Master.Beer  .Horses.Brit     
sol6	h1.Green .Blends     .Coffee.Fish  .Norwegian
	h2.Blue  .Prince     .Water .Cats  .German   
	h3.Yellow.Dunhill    .Milk  .Dogs  .Swede    
	h4.Red   .Blue-Master.Beer  .Horses.Brit     
	h5.White .Pall-Mall  .Tea   .Birds .Dane     
sol7	h1.Green .Pall-Mall  .Coffee.Birds .Norwegian
	h2.Blue  .Prince     .Water .Cats  .German   
	h3.White .Blends     .Milk  .Dogs  .Swede    
	h4.Red   .Blue-Master.Beer  .Horses.Brit     
	h5.Yellow.Dunhill    .Tea   .Fish  .Dane     

$offtext