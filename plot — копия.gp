set terminal png font "Verdana,14" size 1000, 1000
set output "Path.png"
plot '1.txt' w lines lw 2 notitle, '1.txt' pt 7 lc 7 ps 3 notitle
