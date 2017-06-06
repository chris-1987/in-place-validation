rm $1.format
rm $1.format.sa5

./format $1 $1.format

./build $1.format $1.format.sa5

./check $1.format $1.format.sa5
