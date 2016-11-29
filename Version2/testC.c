//test whether C works with R
double myFun(int y, double k){
	int x=0;
	double sum=0;
	while(x<100){
		sum=sum+y*k*x;
		++x;
	}
	return sum;
}

void myFun2(void){

}
