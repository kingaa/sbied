# for help and tutorials see http://rich-iannone.github.io/DiagrammeR/

library('DiagrammeR')

myplot <- "
digraph polio {

  graph [layout=neato, overlap=FALSE] 

  # several 'node' statements
	node [shape = circle,
   		style=filled,
        color = '#B2A389'] // sets as circles
  		Births;

	node [shape = circle,
   		style=filled,
        color = '#F2E3B6'] // sets as circles
  		SB1; SB2; SB3; SB4; SB5; SB6 ; SO;

	node [shape = circle,
   		style=filled,
       color = '#D55D45'] // sets as circles
  		IB; IO;

	node [shape = oval,
   		style=filled,
        color = '#BCDAF2'] // sets as circles
  		Immune;

	node [shape = oval,
   		style=filled,
        color = '#026773', fontcolor=white] // sets as circles
  		Symptomatic;

	node [shape = oval,
   		style=filled,
        color = black, fontcolor=white] // sets as circles
  		Reported;

  # several 'edge' statements
    edge [arrowhead=vee]
	    Births->SB1 SB1->SB2 SB2->SB3 SB3->SB4 SB4->SB5 SB5->SB6 SB6->SO
	    SB1->IB SB2->IB SB3->IB SB4->IB SB5->IB SB6->IB
	    # SO->IO
	    # IO->Immune IB->Immune

    # edge [arrowhead=vee, style=dashed]
	    IO-> Symptomatic Symptomatic->Reported

  # a 'graph' statement
  graph [overlap = true, fontsize = 10, mindist=0.4]
}
"
grViz(myplot)
