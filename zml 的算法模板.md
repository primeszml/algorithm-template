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

#### 常见级数

![image-20240825173137052](C:\Users\zml\AppData\Roaming\Typora\typora-user-images\image-20240825173137052.png)

#### FFT

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
#include <bits/stdc++.h>
using namespace std;
constexpr int P = 998244353;
std::vector<int> rev, roots{0, 1};
int power(int a, int b) {
    int res = 1;
    for (; b; b >>= 1, a = 1ll * a * a % P)
        if (b & 1)
            res = 1ll * res * a % P;
    return res;
}
void dft(std::vector<int> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; ++i) rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
    for (int i = 0; i < n; ++i)if (rev[i] < i)std::swap(a[i], a[rev[i]]);
    if (int(roots.size()) < n) {
        int k = __builtin_ctz(roots.size());
        roots.resize(n);
        while ((1 << k) < n) {
            int e = power(3, (P - 1) >> (k + 1));
            for (int i = 1 << (k - 1); i < (1 << k); ++i)
                roots[2 * i] = roots[i], roots[2 * i + 1] = 1ll * roots[i] * e % P;
            ++k;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; ++j) {
                int u = a[i + j], v = 1ll * a[i + j + k] * roots[k + j] % P, x = u + v;
                if (x >= P) x -= P;
                a[i + j] = x, x = u - v;
                if (x < 0)x += P;
                a[i + j + k] = x;
            }
        }
    }
}
void idft(std::vector<int> &a) {
    int n = a.size();
    std::reverse(a.begin() + 1, a.end());
    dft(a);
    int inv = power(n, P - 2);
    for (int i = 0; i < n; ++i)
        a[i] = 1ll * a[i] * inv % P;
}
struct Poly {
    std::vector<int> a;
    Poly() {}
    Poly(int a0) {if (a0)a = {a0};}
    Poly(const std::vector<int> &a1) : a(a1) {while (!a.empty() && !a.back())a.pop_back();}
    int size() const {return a.size();}
    int operator[](int idx) const {if (idx < 0 || idx >= size())return 0;return a[idx];}
    Poly modxk(int k) const {k = std::min(k, size());return Poly(std::vector<int>(a.begin(), a.begin() + k));}
    Poly divxk(int k) const {
        if (size() <= k)
            return Poly();
        return Poly(std::vector<int>(a.begin() + k, a.end()));
    }
    friend Poly operator+(const Poly a, const Poly &b) {
        std::vector<int> res(std::max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); ++i) {
            res[i] = a[i] + b[i];
            if (res[i] >= P)
                res[i] -= P;
        }
        return Poly(res);
    }
    friend Poly operator-(const Poly a, const Poly &b) {
        std::vector<int> res(std::max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); ++i) {
            res[i] = a[i] - b[i];
            if (res[i] < 0)res[i] += P;
        }
        return Poly(res);
    }
    friend Poly operator*(Poly a, Poly b) {
        int sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot)sz *= 2;
        a.a.resize(sz);
        b.a.resize(sz);
        dft(a.a);
        dft(b.a);
        for (int i = 0; i < sz; ++i)a.a[i] = 1ll * a[i] * b[i] % P;
        idft(a.a);
        return Poly(a.a);
    }
    Poly &operator+=(Poly b) {return (*this) = (*this) + b;}
    Poly &operator-=(Poly b) {return (*this) = (*this) - b;}
    Poly &operator*=(Poly b) {return (*this) = (*this) * b;}
    Poly deriv() const {
        if (a.empty())return Poly();
        std::vector<int> res(size() - 1);
        for (int i = 0; i < size() - 1; ++i)res[i] = 1ll * (i + 1) * a[i + 1] % P;
        return Poly(res);
    }
    Poly integr() const {
        if (a.empty())return Poly();
        std::vector<int> res(size() + 1);
        for (int i = 0; i < size(); ++i)res[i + 1] = 1ll * a[i] * power(i + 1, P - 2) % P;
        return Poly(res);
    }
    Poly inv(int m) const {
        Poly x(power(a[0], P - 2));
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (2 - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    Poly log(int m) const {return (deriv() * inv(m)).integr().modxk(m);}
    Poly exp(int m) const {
        Poly x(1);
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (1 - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    Poly mulT(Poly b) const {
        if (b.size() == 0)return Poly();
        int n = b.size();
        std::reverse(b.a.begin(), b.a.end());
        return ((*this) * b).divxk(n - 1);
    }
    std::vector<int> eval(std::vector<int> x) const {
        if (size() == 0)
            return std::vector<int>(x.size(), 0);
        const int n = std::max(int(x.size()), size());
        std::vector<Poly> q(4 * n);
        std::vector<int> ans(x.size());
        x.resize(n);
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                q[p] = std::vector<int>{1, (P - x[l]) % P};
            } else {
                int m = (l + r) / 2;
                build(2 * p, l, m);
                build(2 * p + 1, m, r);
                q[p] = q[2 * p] * q[2 * p + 1];
            }
        };
        build(1, 0, n);
        std::function<void(int, int, int, const Poly &)> work = [&](int p, int l, int r, const Poly &num) {
            if (r - l == 1) {
                if (l < int(ans.size()))
                    ans[l] = num[0];
            } else {
                int m = (l + r) / 2;
                work(2 * p, l, m, num.mulT(q[2 * p + 1]).modxk(m - l));
                work(2 * p + 1, m, r, num.mulT(q[2 * p]).modxk(r - m));
            }
        };
        work(1, 0, n, mulT(q[1].inv(n)));
        return ans;
    }
};
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n, m;
    std::cin >> n >> m;
    std::vector<int> a(n + 1);
    for (int i = 0; i <= n; ++i)
        std::cin >> a[i];
    std::vector<int> b(m);
    for (int i = 0; i < m; ++i)
        std::cin >> b[i];
    b = Poly(a).eval(b);
    for (int i = 0; i < m; ++i)
        std::cout << b[i] << "\n";
    return 0;
}
```

##### 散装板

```
```



### 生成函数相关



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
    void fomot() {
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

