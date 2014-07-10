mencoder mf://output/*.png -mf w=640:h=480:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
