# zml 的算法模板

## 一、dp模板

### 数位dp

```c++
int dfs(int pos, int pre, int lead, int limit) {
    if (!pos) {
        边界条件
    }
    if (!limit && !lead && dp[pos][pre] != -1) return dp[pos][pre];
    int res = 0, up = limit ? a[pos] : 无限制位;
    for (int i = 0; i <= up; i ++) {
        if (不合法条件) continue;
        res += dfs(pos - 1, 未定参数, lead && !i, limit && i == up);
    }
    return limit ? res : (lead ? res : dp[pos][sum] = res);
}
int cal(int x) {
    memset(dp, -1, sizeof dp);
    len = 0;
    while (x) a[++ len] = x % 进制, x /= 进制;
    return dfs(len, 未定参数, 1, 1);
}
signed main() {
    cin >> l >> r;
    cout << cal(r) - cal(l - 1) << endl;
}
```

## 二、数学模板

### 常见级数

![image-20240825173137052](C:\Users\zml\AppData\Roaming\Typora\typora-user-images\image-20240825173137052.png)

### FFT

```c++
const double pi = acos(-1);
int get_num(int x) // 返回第一个大于等于 x 的 2 的幂。
{
    int res = 1;
    while (res < x) res *= 2;
    return res;
}
void FFT(vector<complex<double>> &a,int n,int op) // 循坏迭代法 FFT
{
    // 蝴蝶变换
    vector<int> re(n,0);
    for (int i = 1; i < n; i++) re[i] = re[i / 2] / 2 + (i & 1) * (n / 2);
    for (int i = 0; i < n; i++) if (i < re[i]) swap(a[re[i]], a[i]);
    // 开始 FFT
    for (int len = 2; len <= n; len *= 2) //, 枚举块长
    {
        for (int i = 0; i < n; i += len) // 枚举块数，块的起点
        {
            complex<double> wk = {1, 0}, w1 = {cos(2 * pi / len), op * sin(2 * pi / len)};
            for (int j = 0; j < len / 2; j++) { // 枚举块内位置
                auto x = a[i + j], y = a[i + j + len / 2] * wk;
                a[i + j] = x + y, a[i + j + len / 2] = x - y;
                wk = wk * w1;
            }
        }
    }
}
struct polynomial_complex {
    int n{}; // n 为 最高次 +1，也就是项数
    vector<complex<double>> a;
    void exlen(int len){ // 扩展长度到 len
        while(a.size() < len) a.push_back({0,0});
    }
    polynomial_complex operator*(polynomial_complex B) {
        polynomial_complex res;
        int sn = get_num(n + B.n - 1);
        res.exlen(sn); exlen(sn); B.exlen(sn);
        FFT(a, sn, 1);
        FFT(B.a, sn, 1);
        for (int i = 0; i < sn; i++) a[i] = a[i] * B.a[i];
        FFT(a, sn, -1);
        res.n = n + B.n - 1;
        for (int i = 0; i < res.n; i++) res.a[i].real((a[i].real() + 0.5) / sn);
        return res;
    }
};
```

### NTT

```c++
void NTT(int a[],int n,int op,int mod){ //n 是项数（需要拓展到 2 的幂），op 为 1 时系数转点值，-1 时点值转系数，mod 是常用模数
    // 蝴蝶变换
    int re[N];
    for(int i=1;i<n;i++) re[i] = re[i/2] / 2 + (i & 1) * (n / 2);
    for(int i=0;i<n;i++) if(re[i] > i) swap(a[i],a[re[i]]);
    // 开始 NTT
    int gi = qmi(g,mod-2,mod);
    for(int len = 2; len <= n; len *= 2){ // 枚举块长
        int g1 = qmi((op == 1? g : gi),(mod-1)/len,mod);
        for(int i = 0; i< n; i += len){ // 枚举块内起点
            int gk = 1;
            for(int j=0;j<len/2;j++){ // 枚举块内位置
                int x = a[i + j], y = (int)a[i + j + len / 2] * gk % mod;
                a[i + j] = (x + y) % mod, a[i + j + len / 2] = ((x + mod - y) % mod);
                gk = (int)gk * g1 % mod;
            }
        }
    }
}
```

### 拓展中国剩余定理

```c++
int exgcd(__int128 a,__int128 b,__int128 &x,__int128 &y){
    if(b == 0){
        x = 1, y = 0;
        return a;
    }
    int d = exgcd(b,a%b,y,x);
    y -= a / b * x;
    return d;
}
bool excrt(__int128 &a1,__int128 &m1,int a2,int m2){ // 拓展中国剩余定理，合并两个同余方程
    __int128 k1,k2;
    __int128 d = exgcd(m1,m2,k1,k2);
    if((a1 - a2) % d) return false; // 表示无解，无法合并
    k1 *= (a1 - a2) / d; k1 %= m2;
    a1 = (a1 - m1 * k1), m1 = m1 / d * m2;
    a1 = (a1 % m1 + m1) % m1;
    return true;
}
```



### 任意模数 NTT

```c++
//
// Created by zml on 2024/8/23.
//
#include<bits/stdc++.h>
#define int long long
using namespace std;

const int N = 3e5 + 10;
const int mod1 = 998244353, mod2 = 167772161, mod3 = 469762049, mod4 = 1004535809; // NTT 常用模数
const int g = 3; // 上面四个模数的原根
int re[N],F1[N],F2[N],F3[N],G1[N],G2[N],G3[N];
int ans[N];
int n,m,p;

int exgcd(__int128 a,__int128 b,__int128 &x,__int128 &y){
    if(b == 0){
        x = 1, y = 0;
        return a;
    }
    int d = exgcd(b,a%b,y,x);
    y -= a / b * x;
    return d;
}

bool excrt(__int128 &a1,__int128 &m1,int a2,int m2){ // 拓展中国剩余定理，合并两个同余方程
    __int128 k1,k2;
    __int128 d = exgcd(m1,m2,k1,k2);
    if((a1 - a2) % d) return false; // 表示无解，无法合并
    k1 *= (a1 - a2) / d; k1 %= m2;
    a1 = (a1 - m1 * k1), m1 = m1 / d * m2;
    a1 = (a1 % m1 + m1) % m1;
    return true;
}

int get_num(int x){
    int temp  = 1;
    while(temp < x) temp *= 2;
    return temp;
}

int qmi(int a,int b,int mod){
    int res = 1;
    a %= mod;
    while(b){
        if(b & 1) res = res * a % mod;
        a = a * a % mod;
        b = b >> 1;
    }
    return res;
}

void NTT(int a[],int n,int op,int mod){
    // 蝴蝶变换
    for(int i=1;i<n;i++) re[i] = re[i/2] / 2 + (i & 1) * (n / 2);
    for(int i=0;i<n;i++) if(re[i] > i) swap(a[i],a[re[i]]);
    // 开始 NTT
    int gi = qmi(g,mod-2,mod);
    for(int len = 2; len <= n; len *= 2){ // 枚举块长
        int g1 = qmi((op == 1? g : gi),(mod-1)/len,mod);
        for(int i = 0; i< n; i += len){ // 枚举块内起点
            int gk = 1;
            for(int j=0;j<len/2;j++){ // 枚举块内位置
                int x = a[i + j], y = (__int128)a[i + j + len / 2] * gk % mod;
                a[i + j] = (x + y) % mod, a[i + j + len / 2] = ((x + mod - y) % mod);
                gk = (__int128)gk * g1 % mod;
            }
        }
    }
}

signed main(){
    cin >> n >> m >> p;
    for(int i=0;i<=n;i++) cin >> F1[i];
    for(int i=0;i<=n;i++) F2[i] = F3[i] = F1[i];
    for(int i=0;i<=m;i++) cin >> G1[i];
    for(int i=0;i<=m;i++) G2[i] = G3[i] = G1[i];
    int sn = n + m + 1;
    sn = get_num(sn);
    int inv1 = qmi(sn,mod1-2,mod1), inv2 = qmi(sn,mod2-2,mod2), inv3 = qmi(sn,mod3 - 2, mod3);
    NTT(F1,sn,1,mod1), NTT(F2,sn,1,mod2), NTT(F3,sn,1,mod3);
    NTT(G1,sn,1,mod1), NTT(G2,sn,1,mod2), NTT(G3,sn,1,mod3);
    for(int i=0;i<sn;i++){
        F1[i] = F1[i] * G1[i] % mod1;
        F2[i] = F2[i] * G2[i] % mod2;
        F3[i] = F3[i] * G3[i] % mod3;
    }
    NTT(F1,sn,-1,mod1), NTT(F2,sn,-1,mod2), NTT(F3,sn,-1,mod3);
    for(int i=0;i<sn;i++) F1[i] = F1[i] * inv1 % mod1, F2[i] = F2[i] * inv2 % mod2, F3[i] = F3[i] * inv3 % mod3;
    for(int i = 0; i < sn ;i ++){
        __int128 a = 0, b = 1;
        excrt(a,b,F1[i],mod1);
        excrt(a,b,F2[i],mod2);
        excrt(a,b,F3[i],mod3);
        ans[i] = a % p;
        if(i <= (n + m))cout << ans[i] << ' ';
    }
    cout << '\n';
}
```

### 多项式求逆

```c++
//
// Created by zml on 2024/8/23.
//
#include<bits/stdc++.h>
#define int long long
using namespace std;
const int N = 4e5 + 10;
const int mod = 998244353, mod2 = 167772161, mod3 = 469762049, mod4 = 1004535809; // NTT 常用模数
const int g = 3; // 上面四个模数的原根
int get_num(int x){
    int temp  = 1;
    while(temp < x) temp *= 2;
    return temp;
}
int qmi(int a,int b){
    int res = 1;
    a %= mod;
    while(b){
        if(b & 1) res = res * a % mod;
        a = a * a % mod;
        b = b >> 1;
    }
    return res;
}
int re[N];
void NTT(int a[],int n,int op){ // n 是项数，将多项式 a 在点值表示和系数表示中转换，点值转系数是注意结果最后还要除以 n，mod 是 NTT 模数
    // 蝴蝶变换
    for(int i=1;i<n;i++) re[i] = re[i/2] / 2 + (i & 1) * (n / 2);
    for(int i=0;i<n;i++) if(re[i] > i) swap(a[i],a[re[i]]);
    // 开始 NTT
    int gi = qmi(g,mod-2);
    for(int len = 2; len <= n; len *= 2){ // 枚举块长
        int g1 = qmi((op == 1? g : gi),(mod-1)/len);
        for(int i = 0; i< n; i += len){ // 枚举块内起点
            int gk = 1;
            for(int j=0;j<len/2;j++){ // 枚举块内位置
                int x = a[i + j], y = a[i + j + len / 2] * gk % mod;
                a[i + j] = (x + y) % mod, a[i + j + len / 2] = ((x + mod - y) % mod);
                gk = gk * g1 % mod;
            }
        }
    }
}
int A[N],B[N];
void mul(int a[],int n,int b[],int m,int c[]){ // a * b 并且将结果存在 c 中，其中 a 的项数是 n, b 的项数是 m，mod 是 模数 (系数的模数)
    int sn = get_num(n + m - 1);
    memset(A,0,sn * 8), memset(B,0,sn * 8);
    for(int i=0;i<n;i++) A[i] = a[i];
    for(int i=0;i<m;i++) B[i] = b[i];
    NTT(A,sn,1), NTT(B,sn,1);
    memset(c,0,sn * 8 ); // 不能提前清空，因为 c 可能和 a、b是同一个
    for(int i=0;i<sn;i++) c[i] = A[i] * B[i] % mod;
    NTT(c,sn,-1);
    int inv = qmi(sn,mod-2);
    for(int i=0;i<sn;i++) c[i] = c[i] * inv % mod;
}
int AA[N],CC[N];
void polyinv(int a[],int n,int b[]) { // 求 a 模 x ^ n 意义下的逆，结果存在 b 中
    int lim = 1;
    memset(AA,0,sizeof(AA)); AA[0] = qmi(a[0],mod-2);
    while(lim < n){
        memset(CC,0,sizeof(CC));
        mul(AA,lim,AA,lim,CC);
        mul(a,n,CC,lim * 2 - 1,CC);
        for(int i=0;i<lim*2;i++) AA[i] = (2 * AA[i] - CC[i] + mod) % mod;
        lim *= 2;
    }
    for(int i=0;i<lim;i++) b[i] = AA[i];
}
int n;
int f[N],ans[N];

signed main(){
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    cin >> n;
    for(int i=0;i<n;i++) cin >> f[i];
    polyinv(f,n,ans);
    for(int i=0;i<n;i++) cout << ans[i] << ' ';
}
```

