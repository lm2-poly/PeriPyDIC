#!/usr/bin/bash
## Tests
cd ../test
echo "--2D direction x+"
python ../pd_dic.py -i input_elas_2D_x+.yaml -t pd > 2D_x+.dat
sed -i '$ d' 2D_x+.res
DIFF=$(diff 2D_x+.dat 2D_x+.res )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi
echo "--2D direction x-"
python ../pd_dic.py -i input_elas_2D_x-.yaml -t pd > 2D_x-.dat
sed -i '$ d' 2D_x-.res
DIFF=$(diff 2D_x-.dat 2D_x-.res )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi
echo "--2D direction y+"
python ../pd_dic.py -i input_elas_2D_y+.yaml -t pd > 2D_y+.dat
sed -i '$ d' 2D_y+.res
DIFF=$(diff 2D_y+.dat 2D_y+.res )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi
echo "--2D direction y-"
python ../pd_dic.py -i input_elas_2D_y-.yaml -t pd > 2D_y-.dat
sed -i '$ d' 2D_y-.res
DIFF=$(diff 2D_y-.dat 2D_y-.res )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi


