#! /bin/bash
cd "/home/marin/Codes/Spin1"
make GA_sp_2.7
for num in {0..10}
do
./GA_sp_2.7  $( bc <<<  $num*0.001*3.14159265359-90*0.001*3.14159265359 ) &

echo $( bc <<<  $num*0.001*3.14159265359-90*0.001*3.14159265359 )
done
exit
