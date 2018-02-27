binwidth = 0.1
set boxwidth binwidth
sum = 0

s(x)          = ((sum=sum+1), 0)
bin(x, width) = width*floor(x/width) + binwidth/2.0

plot "data.txt" u ($1):(s($1))
plot "data.txt" u (bin($1, binwidth)):(1.0/(binwidth*sum)) smooth freq w boxes
#plot "data.txt" u (bin($1, binwidth)):(1.0) smooth freq w boxes
#plot 'data.txt' using 1, '' using 2, '' using 3, '' using 4
set term png 
set output 'plot.png'
replot
set term x11
