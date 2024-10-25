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

### 多项式相关

#### FFT

```c++
#include<bits/stdc++.h>
#define int long long
using namespace std;

const int N = 4e5 + 100;
#define PI acos(-1)
int re[N];
vector<int> up(1,1);
int get_up(int x){
    while(up.size() < x + 1) up.push_back(up.back() >= up.size()? up.back() : up.back() * 2);
    return up[x];
}
void FFT(vector<complex<double>> &a,int op){
    // 蝴蝶变换
    int n = a.size();
    for(int i=1;i<n;i++) re[i] = (re[i>>1] >> 1) + (i & 1) * (n >> 1);
    for(int i=0;i<n;i++) if(i < re[i]) swap(a[re[i]],a[i]);
    // 开始 FFT
    for(int len = 2; len <= n; len *= 2){
        complex<double> w1(cos(2*PI/len),op * sin(2*PI/len));
        for(int i=0;i<n;i+=len){
            complex<double> wk(1,0);
            for(int j=0;j<len/2;j++){
                complex<double> x = a[i + j], y = (wk * a[i + j + len/2]);
                a[i+j] = x + y, a[i + j + len / 2] = x - y;
                wk = wk * w1;
            }
        }
    }
    if(op == -1){
        for(int i=0;i<n;i++) a[i] /= n;
    }
}

struct poly_complex{
    vector<complex<double>> a;
    poly_complex() {}
    poly_complex(complex<double> a1) {
        a = {a1};
    }
    poly_complex(const vector<complex<double>> &a1){a = a1;}
    int size(){return a.size();}
    complex<double> operator[](int idx) const {
        if (idx < 0 || idx >= a.size())return {0,0};
        return a[idx];
    }
    poly_complex modxk(int k) const { // mod x^k
        k = min(k, (int) a.size());
        return poly_complex(vector<complex<double>>(a.begin(), a.begin() + k));
    }
    poly_complex divxk(int k) const { // 除 x ^ k
        if (a.size() <= k)return poly_complex();
        return poly_complex(vector<complex<double>>(a.begin() + k, a.end()));
    }
    friend poly_complex operator+(const poly_complex a, const poly_complex &b) {
        vector<complex<double>> res(max(a.a.size(), b.a.size()));
        for (int i = 0; i < (int)(res.size()); ++i) {
            res[i] = a[i] + b[i];
        }
        return poly_complex(res);
    }
    friend poly_complex operator-(const poly_complex a, const poly_complex &b) {
        vector<complex<double>> res(max(a.a.size(), b.a.size()));
        for (int i = 0; i < (int)(res.size()); ++i) {
            res[i] = a[i] - b[i];
        }
        return poly_complex(res);
    }
    friend poly_complex operator*(poly_complex a, poly_complex b) {
        int sz = a.size() + b.size() - 1;
        sz = get_up(sz);
        a.a.resize(sz);
        b.a.resize(sz);
        FFT(a.a,1), FFT(b.a,1);
        for (int i = 0; i < sz; ++i)a.a[i] = a[i] * b[i];
        FFT(a.a,-1);
        return poly_complex(a.a);
    }
    poly_complex &operator+=(poly_complex b) {return (*this) = (*this) + b;}
    poly_complex &operator-=(poly_complex b) {return (*this) = (*this) - b;}
    poly_complex &operator*=(poly_complex b) {return (*this) = (*this) * b;}
    poly_complex deriv() const { // 求导
        if (a.empty())return poly_complex();
        vector<complex<double>> res(a.size() - 1);
        for (int i = 0; i < a.size() - 1; ++i)res[i] = a[i + 1] * complex<double>(i+1,0);
        return poly_complex(res);
    }
    poly_complex integr() const { // 积分
        if (a.empty())return poly_complex();
        vector<complex<double>> res(a.size() + 1);
        for (int i = 0; i < a.size(); ++i)res[i + 1] = a[i] * complex<double>(1.0/(i+1),0);
        return poly_complex(res);
    }
    poly_complex inv(int m) const { //求逆
        poly_complex x(complex<double>(1,0) / a[0]);
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (poly_complex((2)) - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    poly_complex log(int m) const {return (deriv() * inv(m)).integr().modxk(m);} // 对数
    poly_complex exp(int m) const { // 指数
        poly_complex x(1);
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (poly_complex(1) - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    complex<double> val(double x){
        complex<double> res(0,0), temp(1,0);
        for(auto v : a){
            res = res + v * temp;
            temp = temp * x;
        }
        return res;
    }
    void print(){
        for(auto v : a) cout << v.real() << ' ';
        cout << '\n';
    }
};

void solve(){
    int la,lb,lc;
    double l,r;
    cin >> la >> lb >> lc >> l >> r;
    vector<complex<double>> a(la+1),b(lb+1),c(lc+1);
    for(int i=0;i<=la;i++) {
        double val; cin >> val; a[i] = val;
    }
    for(int i=0;i<=lb;i++) {
        double val; cin >> val; b[i] = val;
    }
    for(int i=0;i<=lc;i++) {
        double val; cin >> val; c[i] = val;
    }
    poly_complex A(a),B(b),C(c);
    A = A * A, B = B * B, C = C * C;
    vector<complex<double>> f(max(A.size(),B.size()),C.size());
    for(int i=0;i<f.size();i++) f[i] = C[i] - A[i] - B[i];
    poly_complex F(f);
    poly_complex G = F.deriv();
    long double x0 = (l + r) / 2;
    while(x0 >= l && x0 <= r){
        long double x = x0 - (F.val(x0) / G.val(x0)).real();
        if(abs(x - x0) < 1e-8) {
            x0 = x;
            break;
        }
        x0 = x;
    }
    if(x0 > l && x0 < r) printf("%.12Lf\n",x0);
    else cout << "Inconsistent!\n";
}

signed main(){
    solve();
}
```

#### NTT

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

#### 拓展中国剩余定理

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



#### 任意模数 NTT

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

#### 多项式全家桶——NTT

##### 高度封装版

```c++
#include<bits/stdc++.h>
#define int long long
using namespace std;

const int N = 3e5 + 100;
const int mod = 998244353, mod2 = 167772161, mod3 = 469762049, mod4 = 1004535809; // NTT 常用模数
const int g = 3; // 上面四个模数的原根

vector<int> rev;
int qmi(int a, int b){
    int res = 1;
    while (b) {
        if (b & 1) res = res * a % mod;
        a = a * a % mod;
        b = b >> 1;
    }
    return res;
}
int gi = qmi(g, mod - 2);
void NTT(vector<int> &a, int op) {
    int n = a.size();
    if (rev.size() != n) {
        rev.resize(n);
        for (int i = 1; i < n; i++) rev[i] = rev[i / 2] / 2 + (i & 1) * (n / 2);
    }
    for (int i = 0; i < n; i++) if (rev[i] < i) swap(a[rev[i]], a[i]);
    for (int len = 2; len <= n; len *= 2) {
        int g1 = qmi((op == 1 ? g : gi), (mod - 1) / len);
        for (int i = 0; i < n; i += len) {
            int gk = 1;
            for (int j = 0; j < len / 2; j++) {
                int x = a[i + j], y = a[i + j + len / 2] * gk % mod;
                a[i + j] = (x + y) % mod, a[i + j + len / 2] = (x + mod - y) % mod;
                gk = gk * g1 % mod;
            }
        }
    }
    if (op == -1) {
        int inv = qmi(n, mod - 2);
        for (int i = 0; i < n; i++) a[i] = a[i] * inv % mod;
    }
}
struct poly {
    vector<int> a;
    poly() {}
    poly(const std::vector<int> &a) : a(a) {} // 传vector构造
    poly(const std::initializer_list<int> &a) : a(a) {}  // 列表构造
    int size() const {
        return a.size();
    }
    void resz(int n) {
        a.resize(n);
    }
    int operator[](int idx) const {
        if (idx < size()) {
            return a[idx];
        } else return 0;
    }
    int &operator[](int idx) {
        return a[idx];
    }
    poly mulxk(int k) { // 返回乘上 x ^ k 的结果，自身值不变
        auto b = a;
        b.insert(b.begin(), k, 0); // 在 b 开头插入 k 个 0 即可
        return poly(b);
    }
    poly modxk(int k) const { // 返回对 x^k 取模的结果，自身值不变
        k = min(k, size());
        return poly(vector<int>(a.begin(), a.begin() + k));
    }
    poly divxk(int k) const { // 返回除以 x ^ k 的结果，自身值不变
        if (k >= size()) {
            return poly();
        }
        return vector<int>(a.begin() + k, a.end());
    }
    friend poly operator+(const poly &a, const poly &b) { // 两个多项式相加，只返回结果，不改变 a,b 的值
        vector<int> res(max(a.size(), b.size()));
        for (int i = 0; i < res.size(); i++) {
            res[i] = (a[i] + b[i]) % mod;
        }
        return poly(res);
    }
    friend poly operator-(const poly &a, const poly &b) { // 两个多项式相-，只返回结果，不改变 a,b 的值
        vector<int> res(max(a.size(), b.size()));
        for (int i = 0; i < res.size(); i++) {
            res[i] = (a[i] - b[i] + mod) % mod;
        }
        return poly(res);
    }
    friend poly operator*(poly a, poly b) { // 返回 a * b 的结果
        if (a.size() == 0 || b.size() == 0) {
            return poly();
        }
        int sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot) sz *= 2;
        a.a.resize(sz), b.a.resize(sz);
        NTT(a.a, 1), NTT(b.a, 1);
        for (int i = 0; i < sz; i++) a[i] = (a[i] * b[i]) % mod;
        NTT(a.a, -1);
        a.a.resize(tot);
        return a;
    }
    friend poly operator*(int a, poly b) {
        for (int i = 0; i < (b.size()); i++) {
            b[i] = a * b[i] % mod;
        }
        return b;
    }
    friend poly operator*(poly a, int b) {
        for (int i = 0; i < (a.size()); i++) {
            a[i] = a[i] * b % mod;
        }
        return a;
    }
    poly &operator+=(poly b) {
        return (*this) = (*this) + b;
    }
    poly &operator-=(poly b) {
        return (*this) = (*this) - b;
    }
    poly &operator*=(poly b) {
        return (*this) = (*this) * b;
    }
    poly &operator*=(int b) {
        return (*this) = (*this) * b;
    }
    poly deriv() const { // 求导
        if (a.empty()) {
            return poly();
        }
        vector<int> res(size() - 1);
        for (int i = 0; i < size() - 1; ++i) {
            res[i] = (i + 1) * a[i + 1] % mod;
        }
        return poly(res);
    }
    poly integr() const { // 不定积分
        std::vector<int> res(size() + 1);
        for (int i = 0; i < size(); ++i) {
            res[i + 1] = a[i] * qmi(i + 1, mod - 2) % mod;
        }
        return poly(res);
    }
    poly inv(int m) const {
        int inv = qmi(a[0], mod - 2);
        poly x{inv};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (poly{2} - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    poly log(int m) const {
        return (deriv() * inv(m)).integr().modxk(m);
    }
    poly exp(int m) const {
        poly x{1};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (poly{1} - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    poly pow(int k,int k1,string str, int m) const {
        // k 是 原始数 % mod, k1 是原始数 % phi（mod），即 mod - 1。str 是原始数
        int i = 0;
        while (i < size() && a[i] == 0) {
            i++;
        }
        if (i == size() || 1LL * i * k >= m) {
            return poly(std::vector<int>(m));
        }
        if(i && (str.size()>7)){
            return poly(std::vector<int>(m));
        }
        int v = a[i];
        auto f = divxk(i) * qmi(v, mod - 2);
        return (f.log(m - i * k) * k).exp(m - i * k).mulxk(i * k) * qmi(v, k1);
    }
    poly sqrt(int m) const {
        poly x{1};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x + (modxk(k) * x.inv(k)).modxk(k)) * ((mod + 1) / 2);
        }
        return x.modxk(m);
    }
};

signed main(){
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    string str;
    int n,k=0,k1=0;
    cin >> n;
    vector<int> a(n);
    for(int i=0;i<n;i++) cin >> a[i];
    poly f(a);
    auto ans = f.inv(n);
    for(int i=0;i<n;i++) cout << ans[i] << ' ';
}
```

##### 散装板（仅整理公式）

```

```

### 线性基

#### 前缀线性基

```c++
struct Prefix_linear_basis {
    int Base[32]{}, pos[32]{};
    int cnt = 0, r{};
    void add(int x, int id) {
        for (int i = 31; i >= 0; i--) {
            if (x >> i & 1) {
                if (Base[i]) {
                    if (pos[i] < id) {
                        swap(Base[i], x);
                        swap(id, pos[i]);
                    }
                    x ^= Base[i];
                } else {
                    Base[i] = x;
                    pos[i] = id;
                    cnt++;
                    break;
                }
            }
        }
    }
    int get_max(int res, int l) {
        for (int i = 31; i >= 0; i--) if (pos[i] >= l) res = max(res, res ^ Base[i]);
        return res;
    }
    int get_min(int res, int l) {
        for (int i = 31; i >= 0; i--) if (pos[i] >= l) res = min(res, res ^ Base[i]);
        return res;
    }
    bool check(int x){
        int res = 0;
        for(int i=31;i>=0;i--){
            if((x >> i & 1) != (res >> i & 1)) res ^= Base[i];
        }
        if(res == x) return true;
        else return false;
    }
    void fomot() { // 求 k-th 前先用这个
        for (int i = 0; i <= 31; i++) {
            for (int j = i + 1; j <= 31; j++) {
                if (Base[j] >> i & 1) Base[j] ^= Base[i];
            }
        }
    }
    int get_kth(int k) {  // 获取第 k 小的数，第 0 小的数记为 0，最小的数能否取 0 取决于基里面的向量个数是否等于 r
        vector<int> tbase;
        fomot();
        for (int i = 0; i < 32; i++)if (Base[i]) tbase.push_back(Base[i]);
        if (tbase.size() == r) k++; // 根据能否取 0 进行调整
        if (k >= (1 << tbase.size())) return -1;
        int res = 0;
        for (int i = 0; i < tbase.size(); i++) {
            if (k >> i & 1) res ^= tbase[i];
        }
        return res;
    }
} p[N];
```

#### 一般线性基

```c++
#include<bits/stdc++.h>
#define int long long
using namespace std;

const int mod = 1e9 + 7;

struct linear_basis {
    int Base[64]{}, cnt, n;
    vector<int> simple_base;

    linear_basis() {
        n = 0, cnt = 0;
        memset(Base, 0, sizeof Base);
        simple_base.clear();
    }

    void add(int x) {
        ++n;
        for (int i = 63; i >= 0; i--) {
            x = min(x, x ^ Base[i]);
        }
        if (!x) return;
        for (int i = 63; i >= 0; i--) {
            Base[i] = min(Base[i], Base[i] ^ x);
        }
        ++cnt;
        for(int i=63;i>=0;i--){
            if(x >> i & 1) {
                Base[i] = x;
                return ;
            }
        }
    }

    int get_max(int res) {
        for (int i = 63; i >= 0; i--) res = max(res, res ^ Base[i]);
        return res;
    }

    int get_min(int res) {
        for (int i = 63; i >= 0; i--) res = min(res, res ^ Base[i]);
        return res;
    }

    int get_kth(int k) {
        simple_base.clear();
        if (n == cnt) k++;
        int res = 0;
        for (int i = 0; i < 64; i++) if (Base[i]) simple_base.push_back(Base[i]);
        if (k >= (1ll << cnt)) return -1;
        for (int i = 0; i < cnt; i++) if (k >> i & 1) res = res ^ simple_base[i];
        return res;
    }

    bool check(int x) {
        for (int i = 63; i >= 0; i--) x = min(x, x ^ Base[i]);
        return x == 0;
    }

    int qmi(int a,int b,int mod){
        int res = 1;
        while(b){
            if(b & 1) res = res * a % mod;
            b = b>> 1;
            a = a * a % mod;
        }
        return res;
    }

    int count(int x) {
        if (!check(x)) return 0;
        return qmi(2, n - cnt, mod);
    }
};

void solve(int test){
    int n;
    cin >> n;
    linear_basis ans;
    for(int i=1;i<=n;i++) {
        int x;
        cin >> x;
        ans.add(x);
    }
    int m;
    cin >> m;
    cout << "Case #" << test << ":\n";
    while(m--){
        int k;
        cin >> k;
        cout << ans.get_kth(k-1) << '\n';
    }
}

signed main(){
    int test;
    cin >> test;
    for(int i=1;i<=test;i++){
        solve(i);
    }
}
```

### 矩阵相关

#### 矩阵快速幂

```c++
const int mod = 998244353;

template<int n, int m>
struct matrix {
    int a[n + 1][m + 1] = {};
    template<int _n>
    matrix operator*(const matrix<m, _n> &B) const {
        matrix<n,_n>res;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= _n; ++j) {
                int sum = 0;
                for (int k = 1; k <= m; ++k) {
                    sum = (sum + a[i][k] * B.a[k][j]) % mod;
                }
                res.a[i][j] = sum;
            }
        }
        return res;
    }
    matrix operator+(const matrix<n, m> &B) const {
        matrix res;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                res.a[i][j] = (a[i][j] + B.a[i][j]) % mod;
            }
        }
        return res;
    }
    void to_one() {
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= m; ++j)
                a[i][j] = (i == j ? 1 : 0);
    }
};
template<int n>
matrix<n, n> qmi(matrix<n, n> A, int b) {
    matrix<n, n> res;
    res.to_one();
    while (b) {
        if (b & 1) res = res * A;
        b = b >> 1;
        A = A * A;
    }
    return res;
}
```

#### 矩阵优化递推常见形式



#### 线段树套矩阵

```c++
//
// Created by zml on 2024/8/2.
//
#include<bits/stdc++.h>
#define int long long
using namespace std;
const int N = 2e5 + 10, mod = 998244353;
template<int n, int m>
struct matrix {
    int a[n + 1][m + 1] = {};
    template<int _n>
    matrix operator*(const matrix<m, _n> &B) const {
        matrix res;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                int sum = 0;
                for (int k = 1; k <= m; ++k) {
                    sum = (sum + a[i][k] * B.a[k][j]) % mod;
                }
                res.a[i][j] = sum;
            }
        }
        return res;
    }
    matrix operator+(const matrix<n, m> &B) const {
        matrix res;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                res.a[i][j] = (a[i][j] + B.a[i][j]) % mod;
            }
        }
        return res;
    }
    void to_one() {
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= m; ++j)
                a[i][j] = (i == j ? 1 : 0);
    }
};
matrix<3, 3> qmi(matrix<3, 3> A, int b) {
    matrix<3, 3> res;
    res.to_one();
    while (b) {
        if (b & 1) res = res * A;
        b = b >> 1;
        A = A * A;
    }
    return res;
}
template<class Node, class Tag, int N>
struct SegmentTree {
#define mid ((l+r)/2)
    Node node[(N + 1) << 2];
    Tag tag[(N + 1) << 2];
    int n;
    void init(int _n) {
        n=_n;
        for (int i = 1; i <= 4 * n; i++)node[i] = Node();
        for (int i = 1; i <= 4 * n; i++)tag[i] = Tag();
    }
    void apply(int p, const Tag &v) {
        tag[p] += v;
        node[p] += v;
    }
    void pushdown(int p) {
        if(tag[p] == Tag()){
            return ;
        }
        apply(2 * p, tag[p]);
        apply(2 * p + 1, tag[p]);
        tag[p] = Tag();
    }
    void pushup(int p) {
        node[p] = node[2 * p] + node[2 * p + 1];
    }
    void updata(int p, int l, int r, int x, int y,const Tag &v) {
        if (l >= x && r <= y) return apply(p,v);
        pushdown(p);
        if (x <= mid)updata(2 * p, l, mid, x, y, v);
        else updata(2 * p + 1, mid + 1, r, x, y, v);
        pushup(p);
    }
    void updata(int l,int r,const Tag&v){
        updata(1,1,n,l,r,v);
    }
    Node query(int p,int l,int r,int ql,int qr){
        if(l >= ql && r <= qr) return node[p];
        if(l > qr || r < ql) return Node();
        pushdown(p);
        return query(p*2,l,mid,ql,qr) + query(p*2+1,mid+1,r,ql,qr);
    }
};
struct Tag {
    matrix<3,3> mul;
    Tag() {mul.to_one();}
    void operator+=(const Tag &v) {
        mul = v.mul;
    }
    bool operator==(const Tag &v){
        for(int i=1;i<=2;i++){
            for(int j=1;j<=2;j++){
                if(mul.a[i][j] != v.mul.a[i][j]) return false;
            }
        }
        return true;
    }
};
struct Node {
    matrix<3,3> sum,suml,sumr,ss;
    Node(){
        
    }
    void operator+=(const Tag &v) {
        sum = v.mul, suml = v.mul, sumr = v.mul, ss = v.mul;
    }
    friend Node operator+(const Node &a, const Node &b) {
        Node res;
        res.sum = a.sum + b.sum;
        res.sum = res.sum + (a.sumr * b.suml);
        res.suml = a.suml + a.ss * b.suml;
        res.sumr = b.sumr + b.ss * a.sumr;
        res.ss = a.ss * b.ss;
        return res;
    }
};
SegmentTree<Node, Tag, N + 1> seg;

signed main() {
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    matrix<3,3> Base,I;
    matrix<3,1> _0;
    _0.a[1][1] = 1, _0.a[2][1] = 1, _0.a[3][1] = 0;
    Base.a[1][1] = 2, Base.a[1][2] = 2, Base.a[1][3] = mod-1;
    Base.a[2][1] = 1, Base.a[2][2] = 0, Base.a[2][3] = 0;
    Base.a[3][1] = 0, Base.a[3][2] = 1, Base.a[3][3] = 0;
    I.to_one();
    int n,m;
    cin >> n >> m;
    seg.init(n);
    for(int i=1;i<=n;i++){
        int x;
        cin >> x;
        Tag temp;
        temp.mul = I + qmi(Base,x);
        seg.updata(i,i,temp);
    }
    while(m--){
        int op,l,r;
        cin >> op >> l >> r;
        if(op == 1){
            Tag temp;
            temp.mul = I + qmi(Base,r);
            seg.updata(l,l,temp);
        }
        else{
            Node it = seg.query(1,1,n,l,r);
            cout << (it.sum.a[3][1] + it.sum.a[3][2]) % mod << '\n';
        }
    }
}


```

### 同余相关

#### exbsgs

```c++
#include<algorithm>
#include<unordered_map>
#include<cmath>
#include<iostream>
#define int long long
using namespace std;
const int inf = -1e9;

int exgcd(int a,int b,int &x,int &y)
{
    if(b==0)
    {
        x=1,y=0;
        return a;
    }
    int d = exgcd(b,a%b,y,x);
    y -= a/b * x;
    return d;
}

int bsgs(int a, int b, int p)
{
    if (1 % p == b % p) return 0;
    int k = (int)sqrt(p) + 1;
    unordered_map<int, int> hash;
    for (int i = 0, j = b % p; i < k; i ++ )
    {
        hash[j] = i;
        j = j * a % p;
    }
    int ak = 1;
    for (int i = 0; i < k; i ++ ) ak =ak * a % p;

    for (int i = 1, j = ak; i <= k; i ++ )
    {
        if (hash.count(j)) return i * k - hash[j];
        j = j * ak % p;
    }
    return inf;
}

int exbsgs(int a,int b,int p)
{
    b = (b % p + p) % p;
    if (1 % p == b % p)return 0;
    int d = __gcd(a, p);
    if (d > 1)
    {
        if (b % d)return inf;
        b /= d, p /= d;
        int x, y;
        int t = exgcd(a / d, p, x, y);
        x = (x % (p / t) + p / t) % (p / t);
        return exbsgs(a, b * x, p) + 1;
    }
    else return bsgs(a, b, p);
}

signed main() {
    int a, b, p;
    while (cin >> a >> p >> b, a || b || p)
    {
        int res = exbsgs(a, b, p);
        if (res < 0)cout << "No Solution\n";
        else cout << res << "\n";
    }
}
```



### 组合相关

#### 卡特兰数——反射容斥

问题：从 $(0,0)$ 出发，每次只能向上或向右走一个单位，且不经过直线 $A:y = x + c_1$ 和 $B:y = x + c_2$，到达 $(n,m)$ 的方案数（$(n,m)$ 在直线 $A$、$B$ 之间）。

不考虑 $A、B$ 的限制，则答案为 $C_{n+m}^n$。

考虑如何减去不合法的方案：不合法的方案必然经过 $A$ 或 $B$，我们按照第一次经过的直线分成 $A$、$B$ 两类。

我们记

- $f(x,y)$：从 $(0,0)$ 出发，第一次经过的直线是 $A$，到达 $(x,y)$ 的方案数。这里的 $(x,y)$ 和 $(0,0)$ 在 $A$ 的两侧。
- $g(x,y)$：从 $(0,0)$ 出发，第一次经过的直线是 $B$，到达$(x,y)$ 的方案数。这里的 $(x,y)$ 和 $(0,0)$ 在 $B$ 的两侧。

任意一条经过 $A$ 的路径，将其首次经过 $A$ 后的部分沿直线 $A$ 轴对称，都是一条从 $(0,0)$ 到达 $(y-c_1,x+c_1)$ 的路径，反过来也成立。

任意一条经过 $B$ 的路径，将其首次经过 $B$ 后的部分沿直线 $B$ 轴对称，都是一条从 $(0,0)$ 到达 $(y-c_2,x+c_2)$ 的路径，反过来也成立。                                                                                                                                                                                                                                                                           

任意一条先经过 $A$ 的路径沿 $A$ 对称后，也先经过 $A$。$B$ 同理。

那么有：

- $f(x,y) = C_{x+y}^{x} - g(y-c_2,x+c_2)$ 。
- $g(x,y) = C_{x+y}^c - f(y-c_1,x+c_1)$。

$f(x,y)$ 和 $g(x,y)$ 可以在 $O(x+y)$ 的时间内递归求出。

例题：https://www.luogu.com.cn/problem/P3266

`ac_code`:

```c++
//
// Created by 35395 on 2024/10/23.
//
#include<bits/stdc++.h>
using namespace std;
#define int long long
const int N = 3e6 + 100, mod = 1e9 + 7;
int f[N], f_inv[N], inv[N];
void init() {
    inv[1] = f[0] = f_inv[0] = 1;
    for (int i = 2; i < N; i++) inv[i] = 1ll * (mod - mod / i) * inv[mod % i] % mod;
    for (int i = 1; i < N; i++) f[i] = 1ll * f[i - 1] * i % mod;
    for (int i = 1; i < N; i++) f_inv[i] = 1ll * f_inv[i - 1] * inv[i] % mod;
}
int C(int a, int b) {
    if (a < 0 || b < 0 || a < b) return 0;
    return 1ll * f[a] * f_inv[b] % mod * f_inv[a - b] % mod;
}
int n, m, c1, c2;
int G(int x, int y);
int F(int x, int y);
signed main() {
    init();
    cin >> n >> m;
    c2 = -(m + 2);
    c1 = 1;
    cout << ((C(n + m + 1 + n, n) - F(n - c1, n + m + 1 + c1) - G(n - c2, n + m + 1 + c2)) % mod + mod) % mod << '\n';
}
int F(int x, int y) {
    if (y < 0 || x < 0) return 0;
    return C(x + y, x) - G(y - c2, x + c2);
}
int G(int x, int y) {
    if (x < 0 || y < 0) return 0;
    return C(x + y, x) - F(y - c1, x + c1);
}
```

