#include <iostream>
#include <ctime>
using namespace std;

int main(){
	int a[30], b[30];
	for(int i=0;i<30;i++){
		a[i] = rand()%10+10;
		b[i] = rand()%10+10;
		cout<<"a: "<<a[i]<<"b: "<<b[i]<<endl;
	}

}