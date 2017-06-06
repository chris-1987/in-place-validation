for i in `seq 100000`
do
	dd if=enwiki_100M of=enwiki_10240 bs=10240 count=1 skip=$i

	./process.sh enwiki_10240

	echo $i
done
