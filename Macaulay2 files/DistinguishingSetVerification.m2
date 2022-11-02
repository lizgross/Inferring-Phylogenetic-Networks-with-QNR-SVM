--The following code verifies that the set of 1126 invariants in S forms a
--permutation invariant distinguishing set for the 24 semi-directed phylogenetic networks.
--The invariants are imported from the text file "Inv1126.txt," but we also show
--how S was generated below.

--The ideals I2, I3, I4, and I5 are the ideals for the tree, 3-cycle, 4-cycle,
--and double-triangle networks numbered 1, 4, 10, and 22. The ideals of the other 
--networks are obtained by permuting these ideals for each different labeling.

restart;

A = (1,1);C = (1,-1);G = (-1,1); T = (-1,-1);

L = 
{(A,A,A,A),(A,A,C,C),(A,C,A,C),(A,C,C,A),
 (A,C,G,T),(C,A,A,C),(C,A,C,A),(C,A,G,T),
 (C,C,A,A),(C,C,C,C),(C,G,A,T),(C,G,C,G),
 (C,G,T,A),(C,C,G,G),(C,G,G,C)};--JC

x = hashTable{A => 0, C => 1, G => 1, T => 1};--JC
y = hashTable{A => 0, C => 1, G => 2, T => 3};


L1 = toList apply(L,i->q_(y#(i#0),y#(i#1),y#(i#2),y#(i#3)));
L2 = {a_0,a_1,b_0,b_1,c_0,c_1,d_0,d_1,e_0,e_1,
      f_0,f_1,g_0,g_1,h_0,h_1,j_0,j_1,k_0,k_1,l_0,l_1};

R = QQ[L1, MonomialOrder=> Weights => {15:1}];

--Inv1126 is the matrix of 1126 invariants
S = value get "Inv1126.txt";
S = toList(flatten entries S);

--Network Ideals for 1, 4, 10, and 22
(I2,I3,I4,I5) = value get "Ideals.txt";


------------------------------------------------------------------------------------------

----Permutations maps of coordinates   ---------------------------------------------------

------------------------------------------------------------------------------------------

  
--(01)

S01 = {q_(0,0,0,0), q_(0,0,1,1), q_(1,0,0,1), q_(1,0,1,0), q_(1,0,2,3),
       q_(0,1,0,1), q_(0,1,1,0), q_(0,1,2,3), q_(1,1,0,0), q_(1,1,1,1),
       q_(1,2,0,3), q_(1,2,2,1), q_(1,2,3,0), q_(1,1,2,2), q_(1,2,1,2)};
  
s01 = map(R,R,matrix{(S01)});



--(02)
S02 = {q_(0,0,0,0), q_(1,0,0,1), q_(0,1,0,1), q_(1,1,0,0), q_(1,2,0,3),
       q_(0,0,1,1), q_(1,0,1,0), q_(1,0,2,3), q_(0,1,1,0), q_(1,1,1,1), 
       q_(0,1,2,3), q_(1,2,1,2), q_(1,2,3,0), q_(1,2,2,1), q_(1,1,2,2)};

s02 = map(R,R,matrix{(S02)});

--(03)

S03 = {q_(0,0,0,0), q_(1,0,1,0), q_(1,1,0,0), q_(0,1,1,0), q_(1,2,3,0),
       q_(1,0,0,1), q_(0,0,1,1), q_(1,0,2,3), q_(0,1,0,1), q_(1,1,1,1),
       q_(1,2,0,3), q_(1,1,2,2), q_(0,1,2,3), q_(1,2,1,2), q_(1,2,2,1)};

s03 = map(R,R,matrix{(S03)});

--(12)
S12 = {q_(0,0,0,0), q_(0,1,0,1), q_(0,0,1,1), q_(0,1,1,0), q_(0,1,2,3),
       q_(1,0,0,1), q_(1,1,0,0), q_(1,2,0,3), q_(1,0,1,0), q_(1,1,1,1), 
       q_(1,0,2,3), q_(1,1,2,2), q_(1,2,3,0), q_(1,2,1,2), q_(1,2,2,1)};

s12 = map(R,R,matrix{(S12)});

--(13)
S13 = {q_(0,0,0,0), q_(0,1,1,0), q_(0,1,0,1), q_(0,0,1,1), q_(0,1,2,3),
       q_(1,1,0,0), q_(1,0,1,0), q_(1,2,3,0), q_(1,0,0,1), q_(1,1,1,1), 
       q_(1,2,0,3), q_(1,2,1,2), q_(1,0,2,3), q_(1,2,2,1), q_(1,1,2,2)};

s13 = map(R,R,matrix{(S13)});

--(23)
S23 =  {q_(0,0,0,0), q_(0,0,1,1), q_(0,1,1,0), q_(0,1,0,1), q_(0,1,2,3),
        q_(1,0,1,0), q_(1,0,0,1), q_(1,0,2,3), q_(1,1,0,0), q_(1,1,1,1), 
        q_(1,2,3,0), q_(1,2,2,1), q_(1,2,0,3), q_(1,1,2,2), q_(1,2,1,2)};

s23 = map(R,R,matrix{(S23)});


------------------------------------------------------------------------------------------

----All Tree Networks    -----------------------------------------------------------------

------------------------------------------------------------------------------------------

J_1 = I2;
J_2 = s12(I2);
J_3 = s02(I2);


------------------------------------------------------------------------------------------

----All 3-cycle Networks    --------------------------------------------------------------

------------------------------------------------------------------------------------------

J_4 = (I3);
J_5 = s12(I3);
J_6 = s13(I3);
J_7 = s02(I3);
J_8 = s03(I3);
J_9 = s02(s13(I3));


------------------------------------------------------------------------------------------

----All 4-cycle Networks    --------------------------------------------------------------

------------------------------------------------------------------------------------------
--1 and 3 at upper left and lower right corners
J_10  = (I4);
J_11  = s13(I4);

--1 and 2 at upper left and lower right corners
J_12  = s23(I4);
J_13  = s12(s23(I4));

--1 and 0 at upper left and lower right corners
J_14  = s03(I4);
J_15  = s01(s03(I4));

-- 0 and 2 at upper left and lower right corners
J_16  = s01(s23(I4));
J_17  = s02(s01(s23(I4)));

-- 0 and 3 at upper left and lower right  corners
J_18  = s01(I4);
J_19 = s03(s01(I4));

-- 2 and 3 at upper left and lower right corners
J_20 = s12(I4);
J_21 = s23(s12(I4));


------------------------------------------------------------------------------------------

----All Double-triangle Networks    ------------------------------------------------------

------------------------------------------------------------------------------------------

J_22  = (I5);
J_23  = s13(I5);
J_24  = s12(I5);


------------------------------------------------------------------------------------------

----Constructing the set S    ------------------------------------------------------------

------------------------------------------------------------------------------------------

PP = set(flatten entries gens J_1);
for i in {4, 10, 22} do { PP = PP +  set(flatten entries gens J_i)}; 

PP = PP + set(flatten entries s01(matrix{toList PP}));
#PP
PP = PP + set(flatten entries s02(matrix{toList PP}));
#PP
PP = PP + set(flatten entries s03(matrix{toList PP}));
#PP
PP = PP + set(flatten entries s12(matrix{toList PP}));
#PP
PP = PP + set(flatten entries s13(matrix{toList PP}));
#PP
PP = PP + set(flatten entries s23(matrix{toList PP}));
#PP
PP = toList PP;

k = 0;
for i from 0 to 1125 do
  for j from 0 to 1125 do 
    if (PP#i == S#j) then k = k + 1
    
--S and PP have the same elements in different order (hence, k = 1126 after running the loop above)
k

------------------------------------------------------------------------------------------

----Verify S is a distinguishing set    --------------------------------------------------

------------------------------------------------------------------------------------------


--prints (i, j, k) where the kth polynomial in S (indexed from 0) is contained in J_i and not in J_j

for i from 1 to 24 do
  for j from 1 to 24 do 
    for k from 0 to 729 do if (S#k % J_i  == 0) and (S#k % J_j != 0) then (print(i, j, k), break)

  
------------------------------------------------------------------------------------------

----Verify S is permutation invariant    -------------------------------------------------

------------------------------------------------------------------------------------------

SS = set(S);
SS = SS + set(flatten entries s01(matrix{toList SS}));
SS = SS + set(flatten entries s02(matrix{toList SS}));
SS = SS + set(flatten entries s03(matrix{toList SS}));
SS = SS + set(flatten entries s12(matrix{toList SS}));
SS = SS + set(flatten entries s13(matrix{toList SS}));
SS = SS + set(flatten entries s23(matrix{toList SS}));

--SS contains all of the 1126 polynomials in S and is invariant under all 
--transpositions.

#SS




