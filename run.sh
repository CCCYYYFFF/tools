data=$1
con=$2
ls $data|while read line
do
	python /share/workspaces/chenyongfeng/m2/Analysis_v2.py -i $data/$line -p $con
done
wait
