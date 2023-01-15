#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef long double ldb;
const int MAXN = 505;
const ldb eps = 1e-13;
const ldb pi = acos(-1.0);

const ldb _h = 1.054571726e-34;
const ldb m_e = 9.10938291e-31;
const ldb _k = 1.0e-31;

ldb mat[MAXN][MAXN];
ldb A[MAXN][MAXN];
int n;
ldb val[MAXN];

struct Poly{
	vector<ldb>f;
	Poly(){f.clear();}
	int size()const {return f.size();}
	void push(ldb x){f.push_back(x);}
	ldb &operator [] (const int &x){
		return f[x];
	}
	ldb operator [] (const int &x)const{
		return f[x];
	}
	friend Poly operator + (const Poly &A,const Poly &B){
		Poly ret;
		int a=A.size(),b=B.size();
		for(int i=0;i<max(a,b);i++)
			ret.push((i<a?A[i]:0)+(i<b?B[i]:0));
		return ret;
	}
	friend Poly operator * (const Poly &A,const ldb &t){
		Poly ret;
		for(int i=0;i<A.size();i++)ret.push(t*A[i]);
		return ret;
	}
	Poly mul(ldb k,ldb b){
		Poly ret;
		int t=size();
		for(int i=0;i<t;i++)
			ret.push(b*f[i]);
		ret.push(0);
		for(int i=0;i<t;i++) ret[i + 1] += k * f[i];
		return ret;
	}
	ldb getans(ldb x){
		ldb ret=0,t = 1.0;
		for(int i=0;i<size();i++,t *= x) ret += f[i] * t;
		return ret;
	}
}f[MAXN];

void solve(){
	for(int i=1;i<n;i++){
		if(A[i+1][i]==0){
			int k=0;
			for(int j=i+1;j<=n&&!k;j++)
				if(A[i+1][j])k=j;
			for(int j=1;j<=n;j++)swap(A[k][j],A[i+1][j]);
			for(int j=1;j<=n;j++)swap(A[j][k],A[j][i+1]);
		}
		ldb qinv = 1.0 / A[i + 1][i];
		for(int j=i+2;j<=n;j++)if(A[j][i] != 0.0){
			ldb d=-1.0*A[j][i]*qinv;
			for(int k=1;k<=n;k++) A[j][k] += d*A[i+1][k];
			d=-d;
			for(int k=1;k<=n;k++) A[k][i+1] += A[k][j]*d;
		}
	}
	f[0].push(1);
	for(int i=1;i<=n;i++){
		ldb t = -1;
		f[i]=f[i-1].mul(-1,A[i][i]);
		for(int j=i-1;j>=1;j--){
			t *= A[j + 1][j];
			f[i]=f[i]+f[j-1]*(t*A[j][i]);
			t *= -1;
		}
	}
}

ldb a[MAXN][MAXN] = {0};

void solve2(){
	for(int i = 1;i <= n;i ++){
		int tmp = i;
		for(int j = i + 1;j <= n;j ++) if(fabs(a[j][i]) > fabs(a[tmp][i])) tmp = j;
		swap(a[i],a[tmp]);
		for(int j = i + 1;j <= n + 1;j ++) a[i][j] = a[i][j] / a[i][i];
		a[i][i] = 1.0;
		for(int j = 1;j <= n;j ++){
			if(i == j) continue;
			for(int k = i + 1;k <= n + 1;k ++){
				a[j][k] -= a[j][i] * a[i][k];
			}
		}
	}
	ldb sum = 0.0;
	for(int i = 1;i <= n;i ++) sum += a[i][n + 1] * a[i][n + 1];
	sum = sqrt(sum);
	for(int i = 1;i <= n;i ++){
		cout << a[i][n + 1] << (i == n ? '\n' : ' ');
	}
}

ldb U[MAXN]; int uid = 0;

void setU(int id){
	if(id == 0){
		for(int i = 1;i <= n;i ++) U[i] = 0.0;
	}
	if(id == 1){
		for(int i = 1;i <= n;i ++) U[i] = 0.5 * _k * (i - (n + 1) / 2) * (i - (n + 1) / 2);
	}
}

int main(){
	cin >> n >> uid;
	setU(uid);
	for(int i = 1;i <= n;i ++){
		for(int j = 1;j <= n;j ++){
			mat[i][j] = 0.0;
		}
	}
	ldb tmp = _h * _h / 2.0 / m_e;
	for(int i = 1;i <= n;i ++){
		mat[i][i - 1] = -1;
		mat[i][i] = 2.0 + U[i] / tmp;
		mat[i][i + 1] = -1;
	}
	for(int i = 1;i <= n;i ++){
		for(int j = 1;j <= n;j ++){
			A[i][j] = mat[i][j] / n / n;
		}
	}
	solve();
	for(int j = 0;j < f[n].f.size();j ++){
		if(j) cout << " + ";
		cout << f[n][j];
		if(j) cout << "x ^ " << j;
	}
	val[1] = (pi * pi * _h * _h / 2.0 / m_e) / tmp;
	for(int i = 1;i <= n;i ++){
		cout << val[i] << " : ";
		for(int j = 1;j <= n;j ++){
			for(int k = 1;k <= n;k ++){
				a[j][k] = -mat[j][k];
				if(j == k) a[j][k] += val[i];
			}
			a[j][n + 1] = 0.0;
		}
		for(int k = 1;k <= n + 1;k ++) a[1][k] = 0.0;
		a[1][1] = 1.0; a[1][n + 1] = 1.0;
		solve2();
	}
	return 0;
}