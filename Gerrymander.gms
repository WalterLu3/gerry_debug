$ontext
input:
state: statename ex. Washington or Wisconsin
num: number of districts for a country
party: party you want to optimize for
optimize option : if you want to optimize Republicans, 1 , otherwise : 0
$offtext
$title gerrymandering model

option mip = osigurobi;
*allow mipstart for osigurobi
option integer4 = 1;
*time limit
option resLim = 10;

set P  parties  / republicans, democrats/;


set nodes 'fips code';

$onExternalInput

*number of districts to be assigned to
scalar DISTRICT_NUM 'district number'/10/;

*if optimizing for republican, the value will be 1, otherwise -1
scalar forRepub /1/;

*fraction of population upper and lower bound for a district
* washington 0.2 / 3.8
* Wisconsin  0.5 / 1.6
scalar lower_fraction 'population lowerbound fraction' /0.5/;
scalar upper_fraction 'population upperbound fraction' /1.6/;
alias (nodes,i,j);

*read in the voting data
parameter num(nodes<,P)/
$ondelim
*$include vote_AZ.csv
*include vote_Washington.csv
$include vote_Wisconsin.csv
$offdelim
/;

* read in the adjacency data
parameter foo(nodes,nodes)/
$ondelim
*$include adjacency_AZ.csv
*include adj_Washington.csv
$include adj_Wisconsin.csv
$offdelim
/;


$offExternalInput

* construct adjacency arc according to adjacency csv
set arcs(nodes,nodes);
arcs(i,j)$(foo(i,j) = 1) = yes;

*declare the district set and setup the district subset according to district_num
set district/d1*d20/;
set sDistrict(district);
sDistrict(district)$(ord(district)<=DISTRICT_NUM) = yes;
alias(sDistrict,d);

* calculate population upper and lower bounds for every district 
Scalar
    Upper
    Lower
;
scalar totalPeople;
totalPeople = sum((I,P), num(I,P)) ;
Lower = lower_fraction * totalPeople / card(d) ;
Upper = upper_fraction * totalPeople / card(d) ;




* some parameter used for statistical use
parameter pop_each_N(nodes);
pop_each_N(nodes) = sum(P,num(nodes,P));
scalar node_num;
node_num = card(nodes);



free variable
    obj;
nonnegative variable
    x(district,i,j) flow from i to j in district d
    pop(district) total population in a given district
    y(district) the product of total population of and whether it wins in district d
;
    
binary variable 
    assign(i,district) 1 if i is assigned to d otherwise 0
    
    root(district,i) 1 if i s the root in d otherwise 0
    assign_1_root_0(i,district) if i is assigned to d and is not a root return 1 otherwise 0
    z(district) 1 if republican win in district d otherwise 0
    both_assign(i,j,district) if both i and j are assigned to a district return 1 otherwise 0
;



* load inital starting point
$ontext
parameter initial_assign(i,district)/
$ondelim
$include initial_Washington.csv
$offdelim
/;
$offtext

*set initial starting value
*assign.fx(i,district) = initial_assign(i,district);


set
    assign_init(i,district) 1 if i is assigned to d otherwise 0(set)
    
    root_init(district,i) 1 if i s the root in d otherwise 0(set)
    assign_1_root_0_init(i,district) if i is assigned to d and is not a root return 1 otherwise 0(set)
    z_init(district) 1 if republican win in district d otherwise 0
    both_assign_init(i,j,district) if both i and j are assigned to a district return 1 otherwise 0
;
$ontext
$GDXIN initial_point.gdx
$load assign_init=assign_init
$load root_init= root_init
$load assign_1_root_0_init= assign_1_root_0_init
$load z_init = z_init
$load both_assign_init = both_assign_init
$GDXIN


*assign.l(i,district) = 1 $assign_init(i,district)
assign.l(i,district) = assign_init(i,district);
root.l(district,i) = root_init(district,i);
assign_1_root_0.l(i,district) = assign_1_root_0_init(i,district);
z.l(district) = z_init(district);
both_assign.l(i,j,district) = both_assign_init(i,j,district);
$offtext
*end initial value setting



$GDXIN initial_point_Wisconsin_1.gdx
$load assign_init
$load root_init
$load assign_1_root_0_init
$load z_init
$load both_assign_init
$GDXIN

*assign.l(assign_init) = 1; (alternative to assignment syntax)

assign.fx(i,district)$assign_init(i,district) = 1;
root.fx(district,i)$root_init(district,i) = 1;
assign_1_root_0.fx(i,district)$assign_1_root_0_init(i,district) = 1;
z.fx(district)$z_init(district) = 1;
both_assign.fx(i,j,district)$both_assign_init(i,j,district) = 1;

scalar original_win;
original_win = sum(district,z.l(district));

equation
*------------------------------connectiviy constraint------------------------------
    objective objective function
    each_assign_one(i) each county is assigned to one district
    one_root_constr(district) there can only be one root for a district
    root_constr(i,district) 
*    root_order_const1(i,j,district) avoid unncessary search for root
*    root_order_const2(i,j,district) avoid unncessary search for root
    avoid_root_search(i,j,district) avoid unncessary search for root
    flow_district_constr_1(district,i,j) only if a county is assigned to a district can it have flow > 0
    flow_district_constr_2(district,i,j) only if a county is assigned to a district can it have flow > 0
*    min_node_district(district)
    assign_1_root_0_constr_1(i,district) logical construction for assign_1_root_0(id)
    assign_1_root_0_constr_2(i,district) logical construction for assign_1_root_0(id)
    root_balance_1(district,i) the network flow for a root node
    root_balance_2(district,i) the network flow for a root node
    non_root_balance_1(district,i) the network flow for a non-root node
    non_root_balance_2(district,i) the network flow for a non-root node
*------------------------------vote calculation constraint------------------------------
    ppl_constr_lo(district) min Number of PPL should be in the interval for every district
    ppl_constr_up(district) max Number of PPL should be in the interval for every district
*    total_population(district)
    decide_win_rep(district) See if in District J republicans are the winner
    decide_win_rep2(district) See if in District J republicans are the winner
*    y_d_const1(district)
*    y_d_const2(district)
*    y_d_const3(district)

;

$ontext
------------------------------connectiviy constraint------------------------------
$offtext
each_assign_one(i)..
    sum(d, assign(i,d)) =e= 1;
    
one_root_constr(d)..
    sum(i, root(d,i)) =e= 1;
    
root_constr(i,d)..
    root(d,i) =l= assign(i,d);
    
*root_order_const1(i,j,d)..
*    assign(i,d) + assign(j,d) - 2*both_assign(i,j,d) =l= 1;

*root_order_const2(i,j,d)..
*    assign(i,d) + assign(j,d) - both_assign(i,j,d) =g= 1;
    
avoid_root_search(i,j,d)$(ord(i) < ord(j))..
    root(d,i) - root(d,j) + 2*both_assign(i,j,d) =l=3;
    
flow_district_constr_1(d,i,j)$(arcs(i,j) or arcs(j,i))..
    x(d,i,j) =l= (node_num)*(assign(i,d));
    
flow_district_constr_2(d,i,j)$(arcs(i,j) or arcs(j,i))..
    x(d,i,j) =l= (node_num)*(assign(j,d));
    
*min_node_district(d)..
*    sum(i,assign(i,d)) =g= low_county_num;

assign_1_root_0_constr_1(i,d)..
    assign(i,d) + (1-root(d,i)) - assign_1_root_0(i,d) =l= 1;
    
assign_1_root_0_constr_2(i,d)..
    assign(i,d) + (1-root(d,i)) - 2*assign_1_root_0(i,d) =g= 0;
    
root_balance_1(d,i)..
    sum(j$(arcs(i,j) or arcs(j,i)), x(d,i,j)) - sum(j$(arcs(i,j) or arcs(j,i)), x(d,j,i)) - (sum(j,assign(j,d))-1) +
    root(d,i)*(node_num) =l= (node_num);
    
root_balance_2(d,i)..
    sum(j$(arcs(i,j) or arcs(j,i)), x(d,i,j) ) - sum(j$(arcs(i,j) or arcs(j,i)), x(d,j,i)) - (sum(j,assign(j,d))-1) +
    root(d,i)*(-node_num) =g= (-node_num);
    
non_root_balance_1(d,i)..
    sum(j$(arcs(i,j) or arcs(j,i)), x(d,i,j) ) - sum(j$(arcs(i,j) or arcs(j,i)), x(d,j,i)) - (-1) +
    (assign_1_root_0(i,d))*(node_num) =l= node_num;
    
non_root_balance_2(d,i)..
    sum(j$(arcs(i,j) or arcs(j,i)), x(d,i,j) ) - sum(j$(arcs(i,j) or arcs(j,i)), x(d,j,i)) - (-1)=g= 0;

$ontext
------------------------------vote calculattion constraint------------------------------
$offtext
ppl_constr_lo(d)..
    sum(i, assign(i,d)*sum(P,num(i,P))) =g= Lower;
    
ppl_constr_up(d)..
    sum(i, assign(i,d)*sum(P,num(i,P))) =l= Upper;

*total_population(d)..
*    sum(i, assign(i,d)*sum(P,num(i,P)))/totalPeople =e= pop(d);
    
decide_win_rep(d)..
    sum(i, assign(i,d)*(num(i,'republicans') - num(i,'democrats'))) - (totalPeople)*z(d) =l= 0;
    
decide_win_rep2(d)..
    sum(i, assign(i,d)*(num(i,'republicans') - num(i,'democrats'))) + (-totalPeople-1)*z(d) =g= -totalPeople ;

*y_d_const1(d)..
*    y(d) =l= z(d);
*y_d_const2(d)..
*    y(d) =l= pop(d);
*y_d_const3(d)..
*    pop(d) - y(d) =l= 1-z(d);

*objective..
*    forRepub * obj =e= 100 * sum(d,y(d));

objective..
    obj =e= forRepub * sum(d,z(d));

model spanning_tree /all/;

*spanning_tree.optfile = 1;

*keep track of solving time
scalar starttime;
starttime = jnow;

solve spanning_tree using mip maximizing obj;

scalar elapsed;
elapsed = (jnow - starttime)*24*3600;


* some parameter used for statistical use
parameter
    afterWinLoss(nodes)
    tempWinLoss(nodes,district)
;

tempWinLoss(i,d)$(assign.l(i,d) + z.l(d) = 2) = 1;
afterWinLoss(i)$(sum(d,tempWinLoss(i,d)) = 1) = 1;


set gerryHeader /result1,result2/;
set gerryHeader2 /result1,RorD/; 
set new_votes_header /num_votes,forRep/;
$onExternalOutput
table assign_result(nodes,district,gerryHeader) 'assigned';
table assign_result2(nodes,district,gerryHeader2) 'assigned';
table assign_result3(nodes,district,gerryHeader2) 'assigned';
parameter total_district_won;
parameter runtime_used;
parameter pop_each_D(district) population;
$offExternalOutput
assign_result(nodes,d,'result1') = assign.l(nodes,d);
assign_result2(nodes,d,'result1') = assign.l(nodes,d);
assign_result2(nodes,d,'RorD')$(assign_result(nodes,d,'result1')=1 and afterWinLoss(nodes)=1) = 1;
assign_result2(nodes,d,'RorD')$(assign_result(nodes,d,'result1')=1 and afterWinLoss(nodes)=0) = -1;
assign_result3(nodes,d,'result1') = assign.l(nodes,d);
assign_result3(nodes,d,'RorD')$(assign_result(nodes,d,'result1')=1 and afterWinLoss(nodes)=1) = 1;
assign_result3(nodes,d,'RorD')$(assign_result(nodes,d,'result1')=1 and afterWinLoss(nodes)=0) = -1;
if (forRepub > 0,
    total_district_won = obj.l;
else
    total_district_won = obj.l + card(d);
);

runtime_used = spanning_tree.resusd;

pop_each_D(d) = sum((i,P),num(i,P)*assign.l(i,d));

$ontext
file da /gerry_mander_%party%_%state%.out/;
da.ap = 1;
put da;
put /'initial_point%num% -------------'/;
put /'original votes:',original_win/;
put /'final votes:',total_district_won/;
put /'best possible solution:',spanning_tree.objest/;
put /'solving time:',elapsed/;
putclose da;
$offtext


*display original_win,spanning_tree.solvestat,spanning_tree.rescalc,spanning_tree.sysvers;