#include <bits/stdc++.h>

using namespace std;

#define int long long
#define pb push_back
#define ff first
#define ss second

#define FOR(i, j, k, in) for (int i=j ; i<k ; i+=in)
#define RFOR(i, j, k, in) for (int i=j ; i>=k ; i-=in)
#define REP(i, j) FOR(i, 0, j, 1)
#define RREP(i, j) RFOR(i, j, 0, 1)
#define ALL(v) v.begin(), v.end()
#define SORT(v) sort(ALL(v))
#define REVERSE(v) reverse(ALL(v))
#define SZ(v) (int)v.size()

typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<string> vs;
typedef vector<pii> vii;
typedef vector<vi> vvi;
typedef map<int,int> mpii;

void solve(){
}
signed main() 
{
  int t=1;
  // cin>>t;
  
  while(t--){
    solve();
  }
  return 0;
  
}

/*


Start using personal laptop with different browsers while using different codeforces accounts 


*/
/*https://onecompiler.com/cpp/43u9qm9mz*/
/*https://codeforces.com/problemset/status/17/problem/A*/
/*https://codeforces.com/problemset/problem/116/B*/


/* Codeforces Handle - randomGuyTraining
(try to solve without peaking at the solution) 
(further use 2 laptops to see the solution if no ided what to do)
(login into the second laptop using handle monkeyDboa)*/
/*Alternatively if the question is done in python do it in cpp using same handle*/
/* Codeforces Handle - monkeyDboa */

/*                      */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/* PASTE TILL HERE ONLY */
/*                      */

/* firstever dp problem*/
// const int maxn = 210;
// bool dp[maxn];
// int q,d;
// void solve(){
//   memset(dp,0,sizeof dp);
//   dp[0] = 1;
//   cin>>q>>d;
//   if(!d) d+=10;
//   int mx = d*10;
//   for(int i=0;10*i+d<=mx;i++)
//     for(int j=0;10*i+d+j<=mx;j++)
//       dp[10*i+d+j]|=dp[j];
//   while(q--){
//     int u;
//     cin>>u;
//     if(u>=mx or dp[u])cout<<"YES\n";
//     else cout<<"NO\n";
//   }
// }
 
/*equal prefix sums from 2 arrays*/
// int n,m,c,i,a[100005],b[100005],p[100005];
// void solve(){
//   int n,m;cin>>n>>m;
//   REP(i,n)cin>>a[i+1];
//   REP(i,m)cin>>b[i+1];
//   for(i =2;i<=n;i++)a[i]+=a[i-1];
//   for(i =2;i<=m;i++)b[i]+=b[i-1];
//   for(i =1;i<=n;i++)p[a[i]]++;
//   for(i =1;i<=m;i++)c+=p[b[i]];
//   cout<<c<<endl;
// }

/* function to rotate matrix 90 deg Clockwise*/
// vector<vector<int>> rotate90Clockwise(const vector<vector<int>>& mat) {
//     int n = mat.size();
//     int m = mat[0].size();
//     vector<vector<int>> rotated(m, vector<int>(n));

//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < m; j++) {
//             rotated[j][n - i - 1] = mat[i][j];
//         }
//     }
//     return rotated;
// }
/*function to find mindigit and maxDig*/
// int minDig(int x){
//   int ans = INT_MAX;
//   while(x){
//     ans = min(ans,x%10);
//     x/=10;
//   }
//   return ans;
// }

// int maxDig(int x){
//   int ans = INT_MIN;
//   while(x){
//     ans = max(ans,x%10);
//     x/=10;
//   }
//   return ans;
// }
/*follow string of directions to reach from start to end*/
// char arr[60][60];
// int f[4]={0,1,2,3};
// void solve(){
//   int n,m,x,y;
// 	cin>>n>>m;
// 	for(int i=1;i<=n;i++){
// 		for(int j=1;j<=m;j++){
// 			cin>>arr[i][j];
// 			if(arr[i][j]=='S'){
// 				x=j;y=i;
// 			}
// 		}
// 	}
// 	string str;
// 	cin>>str;
// 	int len = str.length();
// 	int cnt = 0;
// 	do{
// 	  char a = f[0]+'0',b=f[1]+'0',c=f[2]+'0',d=f[3]+'0';
// 	  int sx=x,sy=y;
// 	  REP(i,len){
// 	    if(str[i]==a)sx--;
// 	    if(str[i]==b)sx++;
// 	    if(str[i]==c)sy--;
// 	    if(str[i]==d)sy++;
// 	    if(arr[sy][sx]=='E'){
// 	      cnt++;
// 	      break;
// 	    }
// 	    else if (arr[sy][sx]=='#' or sx>m or sx<1 or sy<1 or sy>n)break;
	    
// 	  }
// 	}while(next_permutation(f,f+4));
//   cout<<cnt<<endl;
// }

/*minimum number of operations from a to x till n is min(n,a-x)*/
// void solve() {
//   int n,a,b,x,y;
//   cin>>a>>b>>x>>y>>n;
//   int ans = 1e18;
//   REP(i,2){
//     int da = min(n,a-x);
//     int db = min(n-da,b-y);
//     ans = min(ans,(a-da)*(b-db));
//     swap(a,b);
//     swap(x,y);
//   }
//   cout<<ans<<endl;
// }
/* function to check the condition between neighbours*/
// int check(vi a){
//   int n = SZ(a);
//   REP(i,n){if(a[i]>=a[(i+1)%n]+a[(i-1+n)%n])return 0;}
//   return 1;
// }
/* stack pointer inc and dec*/
// void solve() {
//   string s;cin>>s;
//   int stk = 0;
//   REP(i,SZ(s)){
//     if(stk and s[i]=='B')stk--;
//     else stk++;
//   }
//   cout<<stk<<endl;
// }
/* check the first occurrence of each alphabet is in sorted order*/
// void solve() {

//   string s;
//   cin>>s;
//   map<char,int> was;
//   for(int i=SZ(s)-1;i>=0;i--)
//     was[s[i]] = i;
//   vector<pair<char,int>> arr;
//   for(auto [k,v]:was)arr.pb({k,v});
//   SORT(arr);
//   vi barr;
//   REP(i,SZ(arr))barr.pb(arr[i].ss);
//   vi carr;
//   REP(i,SZ(arr))carr.pb(arr[i].ff);
//   REP(i,SZ(carr)){if(carr[i]-'a'!=i){cout<<"NO\n";return;}}
//   if(is_sorted(barr.begin(),barr.end()))cout<<"YES\n";
//   else cout<<"NO\n";
// }

/*increment max no from remaining and position of max_element*/
// void solve(){
//   int n;
//   cin>>n;
//   vi v(n);
//   REP(i,n)cin>>v[i];
//   int cnt = 0;
//   while(1){
//     int pos = max_element(v.begin()+1,v.end())-v.begin();
//     if(v[0]>v[pos])break;
//     v[0]++;
//     v[pos]--;
//     cnt++;
    
//   }
//   cout<<cnt<<endl;
// }


/*use of multiset properly*/
// void solve(){
//   int n;
//   cin>>n;
//   multiset<int,greater<int>> s;
//   while(n--){
//     int d;
//     cin>>d;
//     s.insert(d);
//   }
//   int k =*s.begin();
//   cout<<k<<" ";
//   for(int i=1;i<=k;i+=1)if(k%i==0)s.erase(s.find(i));
//   cout<<*s.begin();
// }

/*Lexicographically max subsequence*/
// void solve(){
//   string a;cin>>a;
//   int n = a.size();
//   string ans;
//   REP(i,n){
//     while(ans.size() and a[i]>ans.back())ans.pop_back();
//     ans+=a[i];
//   }
//   cout<<ans<<endl;
// }

/*how to initialize the variable*/
// void solve(){
// int n,x;
// cin>>n>>x;
// int t(1);
// int ans(0);
// REP(i,n){
//   int l,r;
//   cin>>l>>r;
//   t += x*((l-t)/x);
//   ans += r-t+1;
//   t = r+1;
// }
// cout<<ans<<endl;
// }

/*function to check if number is divisible by its digits*/
// bool isfair(int n){
//   int ts = n;
//   while(ts>0){
//     if((ts%10)!=0){
//       if(n%(ts%10))return false;
//     }
//     ts = ts/10;
//   }
//   return true;
// }

/*loop to find all substrings of a string*/
// REP(i,n){
//     FOR(j,i,n,1){
//       string tmp = s.substr(i,j-i+1);
//       cout<<tmp<<endl;
//       
//     }
//   }
  
/*function to check if robot sequence return to start*/
// int check(string s){
//   int cnt1=0,cnt2=0;
//   REP(i,s.size()){
//     cnt1+=s[i]=='U';
//     cnt1-=s[i]=='D';
    
//     cnt2+=s[i]=='L';
//     cnt2-=s[i]=='R';
//   }
//   return cnt1==0 and cnt2==0;
// }

/*prefix and suffix sums*/

// void solve(){
//   int n;
//   cin>>n;
//   string s;
//   cin>>s;
//   vi pref(n+1,0);
//   vi suff(n+1,0);
//   for(int i=1;i<=n;i++){  
//     pref[i] = pref[i-1]+(s[i-1]=='A' or s[i-1]=='F');
//     // cout<<pref[i]<<" ";
//   }
//   // cout<<endl;
//   for(int i = n-2;i>=0;i--){
//     suff[i] = suff[i+1]+(s[i+1]=='A' or s[i+1]=='F');
//     // cout<<suff[i]<<" ";
//   }
//   // cout<<endl;
//   int ans = 0;
//   for(int i=0;i<n;i++){
//     if(s[i]!='F'){
//       int x=0;
//       x=pref[i];
//       int y= suff[i];
//       // cout<<"i="<<i<<",x="<<x<<",y="<<y<<endl;
//       if(x+y==n-1)ans++;
//     }
//   }
//   cout<<ans<<endl;
  
// }
/*convert into lower case*/
// std::transform(str.begin(), str.end(), str.begin(),
//   [](unsigned char c){ return std::tolower(c); });

/*function to convert integer into words*/
// string convertToWords(int n) {
//     if (n == 0) 
//         return "Zero";
    
//     // Words for numbers 0 to 19
//     vector<string> units = {
//         "",        "One",       "Two",      "Three",
//         "Four",    "Five",      "Six",      "Seven",
//         "Eight",   "Nine",      "Ten",      "Eleven",
//         "Twelve",  "Thirteen",  "Fourteen", "Fifteen",
//         "Sixteen", "Seventeen", "Eighteen", "Nineteen"
//     };
    
//     // Words for numbers multiple of 10        
//     vector<string> tens = { 
//         "",     "",     "Twenty",  "Thirty", "Forty",
//         "Fifty", "Sixty", "Seventy", "Eighty", "Ninety" 
//     };
    
//     vector<string> multiplier = 
//     				{"", "Thousand", "Million", "Billion"};
  
//     string res = "";
//     int group = 0;
    
//     // Process number in group of 1000s
//     while (n > 0) {
//         if (n % 1000 != 0) {
            
//             int value = n % 1000;
//             string temp = "";
            
//             // Handle 3 digit number
//             if (value >= 100) {
//                 temp = units[value / 100] + " Hundred ";
//                 value %= 100;
//             }

//             // Handle 2 digit number
//             if (value >= 20) {
//                 temp += tens[value / 10] + " ";
//                 value %= 10;
//             }

//             // Handle unit number
//             if (value > 0) {
//                 temp += units[value] + " ";
//             }

//             // Add the multiplier according to the group
//             temp += multiplier[group] + " ";
            
//             // Add this group result to overall result
//             res = temp + res;
//         }
//         n /= 1000;
//         group++;
//     }
    
//     // Remove trailing space
//     return res.substr(0, res.find_last_not_of(" ") + 1);
// }

/* maximum number of a digit taken from available count*/
// void solve(){
//   int n;
//   cin>>n;
//   int r = 0, t= 0;
//   REP(i,n){
//     int x;
//     cin>>x;
//     if(x==5)r++;
//     else t++;
//   }
//   if(t==0)cout<<"-1\n";
//   else if(r<9){
//     cout<<"0\n";
//   }else{
//     r-=r%9;
//     REP(i,r)cout<<5;
//     REP(i,t)cout<<0;
//     cout<<endl;
//   }
// }

/*function to get all lucky numbers made of 4 and 7*/
// const int N = 1000;
// set<int> luckyNo;
// void getallLuckyNumbers(int x){
//   if(x>N)return;
//   luckyNo.insert(x);
  
//   getallLuckyNumbers(x*10+4);
//   getallLuckyNumbers(x*10+7);
// }

/*function to find sum of digits*/
// int sumOfDigits(int x){
//   int rv = 0;
//   while(x){
//     rv+=x%10;
//     x/=10;
//   }
//   return rv;
// }

/*getting primes using sieve*/
// vector<bool> prime(100010, true);
// vector<int> res;
// void sieve() {
//     int n = 100000;
//     // creation of boolean array
    
//     for (int p = 2; p * p <= n; p++) {
//         if (prime[p] == true) {
            
//             // marking as false
//             for (int i = p * p; i <= n; i += p)
//                 prime[i] = false;
//         }
//     }
    
    
//     for (int p = 2; p <= n; p++){
//         if (prime[p]){ 
//             res.push_back(p);
//         }
//     }
    
// }
/*first put a then put b then put c in a loop*/
// void solve(){
//   int n,k;
//   cin>>n>>k;
//   string s;
//   cin>>s;int div1=n/k;
//   vector<int>freq(26,0);
//   for(int i=0;i<s.length();i++){
//       freq[s[i]-'a']++;
//   }
//   while(k--){
//     char ch;
//     // first put a then b then c and continue till freq[i]>0
//     for(  ch='a';ch<='a'+div1-1;ch++){
//         int i=ch-'a';
//         if(freq[i]>0){
//             freq[i]--;
//         }else{
//             break;
//         }
//     }
//     cout<<ch;
//   }
//   cout<<endl;

// }

/*application of double in division*/
// void solve(){
//   double n,k;
//   cin>>n>>k;
//   vector<double> v(n);
//   double sum=0;
//   for(auto &it: v){ cin>> it; sum += it; }
//   int cnt = n;
//   while(sum/cnt < k-0.5)sum+=k,cnt++;
//   cout<<cnt-n<<endl;
// }

/* finding min element in vector*/
// void solve(){
//   int n,m;
//   cin>>n>>m;
//   vi cnt(n);
//   REP(i,m){
//     int col;
//     cin>>col;
//     cnt[col-1]++;
//   }
//   cout << *min_element(ALL(cnt)) << endl;
  
// }

/*application of back and pop_back in string*/
// void solve(){
//   int n;cin>>n;string s;cin>>s;
//   string res = "";
//   while(!s.empty()){
//     int x;
//     if(s.back()=='a' or s.back()=='e')x=2;
//     else x=3;
//     while(x--){
//       res+=s.back();
//       s.pop_back();
//     }
//     res+='.';
//   }
//   res.pop_back();
//   REVERSE(res);
//   cout<<res<<endl;
// }
/*2d count hashmap related code*/
// int n,m,i,j,a,r,h[2000][5];
// char c;
// void solve(){
//   for(cin>>n>>m;i<n;i++)
//     for(j=0;j<m;j++)
//       cin>>c,h[j][c-'A']++;
// 	for(i=0;i<m;i++)
// 	  cin>>a,r+=a**max_element(h[i],h[i]+5);
// 	 cout<<r<<endl;
// }

/* function to check if a character is vowel or not*/
// bool isVowel(char c){
//   for(char v:"aouie")if(c==v)return 1;
//   return 0;
// }
/* s[i]==s[i+1] means crossing a dot with equation y=x RR or UU or LL or DD*/

/*function to check if perfect square or not*/
// bool isNotPerfectSquare(int n){
//   int y = sqrt(n);
//   return y*y!=n;
// }

/*function to remove trailing zeros in a integer string*/
// string removeTrailingZeros(string s){
//   string t;
//   reverse(s.begin(),s.end());
//   int i=0;
//   while(s[i]=='0')i++;
//   for(int k=i;k<s.size();k++)t+=s[k];
//   reverse(t.begin(),t.end());
//   return t;
// }

/*function to find one of clockwise or ccw direction
 if any one is the answer
*/
// void solve(){
//   char s,e;
//   int t;
//   cin>>s>>e>>t;
//   map<char,int> was = {{'^', 0}, {'>', 1}, {'v', 2}, {'<', 3}};
//   if((was[s]-was[e])%2==0)
//     cout<<"undefined";
//   else if((t%4 + was[s])%4 == was[e])
//     cout<<"cw";
//   else cout<<"ccw";
// }
/*function for cyclicShift of a string*/
// string cyclicShift(string s){
//   int n = s.size();
//   string firstPart = s.substr(0,n-1);
//   string secondPart = s.substr(n-1);
//   return secondPart+firstPart;
// }
/* how to erase a part of string and replace a part of string */
// void solve(){
//   int n;
//   cin>>n;
//   string s;
//   cin>>s;
//   while(s.find("ogo")!=string::npos){
//     int i = s.find("ogo");
//     while(s[i+3]=='g' and s[i+4]=='o')
//       s.erase(i+2,2);
//     s.replace(i,3,"***");
//   }
//   cout<<s;
// }
/*function to remove a element in vector*/
// vi removeEle(vi &a,int idx){
//   vi b;
//   REP(i,a.size())if(i!=idx)b.pb(a[i]);
//   return b;
// }

/*function to remove leading zeros in a integer
  represented as a string
*/
// string removeL(string s){
//   int i=0;
//   while(s[i]=='0'){
//     s[i]='?';
//     i++;
//   }
//   string ss;
//   for(int i=0;i<s.size();i++)if(s[i]!='?')ss+=s[i];
//   return ss;
// }

/* to calculate operations in while(x<d)x+=t efficiently*/
// void solve(){
//   int n,d;
//   cin>>n>>d;
//   vi b(n);
//   REP(i,n)cin>>b[i];
//   int op = 0;
//   for(int i=0;i+1<n;i++){
//     if(b[i+1]<=b[i]){
//       op += (b[i]-b[i+1])/d + 1;
//       b[i+1] += ((b[i]-b[i+1])/d + 1)*d;
//     }
//     // while(b[i+1]<=b[i])b[i+1]+=d,op++;
//   }
//   cout<<op<<endl;
  
// }

/*upper_bound and lower_bound on an array*/
// int a[300005];
// void solve(){
//   int x,y;
//   int n;
//   cin>>n>>x>>y;
//   int sum=0;
//   for(int i=1;i<=n;i++)
//     cin>>a[i],sum+=a[i];
//   sort(a+1,a+1+n);
//   int ans=0;
//   for(int i=1;i<=n-1;i++){
//     // upper bound is taken to exclude one more element
//     ans+=(upper_bound(a+i+1,a+1+n,sum-a[i]-x)-lower_bound(a+i+1,a+1+n,sum-a[i]-y));
//   }
//   cout<<ans<<'\n';
// }

/*Two pointer + prefix sum*/
// void solve(){
//   int n;cin>>n;
//   int a[n+1],h[n+1]; h[0] = 0;
//   FOR(i,1,n+1,1){
//     cin>>a[i];
//     h[i]=h[i-1]+a[i];
//   }
//   string s;
//   cin>>s;
//   int ans = 0;
//   int l = 0, r = n-1;
//   while(l<r){
//     while(l<r and s[l]=='R')l++;
//     while(l<r and s[r]=='L')r--;
//     if(l<r)ans += h[r+1]-h[l];
//     l++,r--;
    
//   }
//   cout<<ans<<endl;
// }

/*Two pointer without overlapping i and j*/
// void solve(){
//   int n,l,r;
//   cin>>n>>l>>r;
//   int a[n];
//   REP(i,n)cin>>a[i];
//   int ans = 0;
//   int sum = 0;
//   for(int i=0,j=0;i<n;i++){
//     sum+=a[i];
//     while(sum>r)sum-=a[j++];
//     if(sum>=l and sum<=r)ans++,j=i+1,sum=0;
//   }
//   cout<<ans<<endl;
// }

/*Two pointer in forward direction*/
// void solve(){
//   int n,s;
//   cin>>n>>s;
//   int len = -1,a[n],sum = 0;
//   for(int i=0,j=0;i<n;i++){
//     cin>>a[i];
//     sum+=a[i];
//     if(sum>s)sum-=a[j++];
//     if(sum==s)len = max(len,i-j+1);
//   }
//   int ans = n-len;
//   if(sum>len or len<0)ans = -1;
//   cout<<ans<<endl;
// }

/*Transition of states using queue*/
// const int N = 5e3 + 50;
// int a[N],vis[N];
// int n;
// void solve(){
//   cin>>n;
//   vis[1] = 1;
//   FOR(i,1,n+1,1)cin>>a[i];
//   vii v;
//   queue<int> q;
//   if(a[1]!=0){
//     q.push(1);
//     vis[1] = 1;
//   }
//   while(q.size()){
//     int x = q.front();
//     int mx = -1;
//     int curr = -1;
//     FOR(i,1,n+1,1){
//       if(vis[i]==1)continue;
//       if(a[i]>mx){
//         mx = a[i];
//         curr = i;
//       }
//     }
//     if(mx==-1)break;
//     vis[curr] = 1;
//     if(a[curr]!=0)q.push(curr);
//     a[x]--;
//     if(a[x]==0)q.pop();
//     v.pb({x,curr});
//   }
//   FOR(i,1,n+1,1){
//     if(vis[i]==0){cout<<"-1\n";return;}
//   }
//   cout<<v.size()<<endl;
//   for(auto [l,r]:v)cout<<l<<" "<<r<<endl;
// }

/* two pointer to iterate through equal elements and minimize the cost */
// void solve(){
//   int n,m,x,y;
//   cin>>n>>m>>x>>y;
//   int ans = 0;
//   y =min(y,2*x);
//   REP(I,n){
//     string s;cin>>s;
//     int i=0;
//     while(i<m){
//       if(s[i]=='*'){
//         i++;
//         continue;
//       }
//       int j=i;
//       while(j+1<m and s[j+1]=='.')j++;
//       int l = j-i+1;
//       ans += l%2 * x + l/2 * y;
//       i= j+1;
//     }
//   }
//   cout<<ans<<endl;
  
// }

/*BFS find any last visited coordinate using bfs*/
// int n, m, k;
// queue <pii> t;
// bool vis[3005][3005];
// pii last;
  
// void solve(){
//   cin >> n >> m >> k;
// 	for (int i = 0; i < k; i++)
// 	{
// 		cin >> last.first >> last.second;
// 		t.push(last);
// 	}
// 	while (!t.empty()){
// 		int x = t.front().first, y = t.front().second;
// 		t.pop();
// 		if (x <= 0 or y <= 0 or x > n or y > m or vis[x][y]) continue;
// 		vis[x][y] = 1;
// 		last = {x, y};
// 		int d[4][2] = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};
// 		for (int i = 0; i < 4; i++){
// 			t.push(make_pair(x + d[i][0], y + d[i][1]));
// 		}
// 	}
// 	cout << last.first << ' ' << last.second << endl;
// }

/*use of memset*/
// void solve(){
//   // cout<<N<<endl;
//   int n;cin>>n;int p[n+1];
//   FOR(i,1,n+1,1)cin>>p[i];
//   int used[n+1];
//   memset(used,0,sizeof used);
//   int rv=0;
//   FOR(i,1,n+1,1){
//     if(!used[i]){
//       int curr = i;
//       int le = 0;
//       while(used[curr]==0){
//         le++;
//         used[curr]=1;
//         curr= p[curr];
//       }
//       rv += (le-1)/2;
//     }
//   }
//   cout<<rv<<endl;
// }
/*find Data structure from UFDS*/
// int mod = (int)1e9 + 7;
// int p[1000010];
// int find(int x){
//   if(p[x]!=x)p[x] = find(p[x]);
//   return p[x];
// }
// int a[1000010],b[1000010],c[1000010],f[1000010];
// void solve(){
//   int n;cin>>n;
//   for(int i=1;i<=n;i++){
//     cin>>a[i];
//     // represents visited array
//     f[i] = 1;
//     // represents parent of each set
//     p[i] = i;
//   }
//   for(int i=1;i<=n;i++){
//     cin>>b[i];
//     p[find(a[i])] = find(b[i]);
//   }
//   for(int i=1;i<=n;i++){
//     cin>>c[i];
//     if(c[i] or a[i]==b[i])f[find(a[i])]=0;
//   }
//   int ans = 1;
//   for(int i=1;i<=n;i++){
//     if(f[i] and find(i)==i)ans = ans*2%mod;
//   }
//   cout<<ans<<endl;
// }

/*BFS marking components and calculating distance using bfs in a 2d grid */
// void solve(){
//   int n;
// 	cin >> n;
// 	int x1, y1, x2, y2;
// 	cin >> x1 >> y1 >> x2 >> y2;
// 	vector<vector<char>>s(n + 1, vector<char>(n + 1));
// 	for(int i = 1; i <= n; i++)
// 	{
// 		for(int j = 1; j <= n; j++)
// 		{
// 			cin >> s[i][j];
// 		}
// 	}
// 	vector<int>dx{-1,1,0,0};
// 	vector<int>dy{0,0,-1,1};
// 	auto check = [&](int x, int y)
// 	{
// 		return x >= 1 && x <= n && y >= 1 && y <= n;
// 	};
// 	vector<vector<int>>dis1(n + 1, vector<int>(n + 1, -1));
// 	auto dis2 = dis1;
// 	auto bfs = [&](int x, int y, vector<vector<int>>&dis)->void
// 	{
// 		queue<pair<int,int>>q;
// 		dis[x][y] = 0;
// 		q.push({x, y});
// 		while(!q.empty())
// 		{
// 			auto [u, v] = q.front();
// 			q.pop();
// 			for(int i = 0; i < 4; i++)
// 			{
// 				int xx = u + dx[i], yy = v + dy[i];
// 				if(check(xx, yy) && s[xx][yy] == '0' && dis[xx][yy] == -1)
// 				{
// 					dis[xx][yy] = dis[u][v] + 1;
// 					q.push({xx, yy});
// 				}
// 			}
// 		}
// 	};
// 	bfs(x1, y1, dis1);
// 	bfs(x2, y2, dis2);
// 	int ans = 1e9;
// 	for(int i = 1; i <= n; i++)
// 	{
// 		for(int j = 1; j <= n; j++)
// 		{
// 			for(int k = 1; k <= n; k++)
// 			{
// 				for(int z = 1; z <= n; z++)
// 				{	
// 					if(dis1[i][j] != -1 && dis2[k][z] != -1)
// 					{
// 						ans = min(ans, (i - k) * (i - k) + (j - z) * (j - z));
// 					}
// 				}
// 			}
// 		}
// 	}
// 	cout << ans << "\n";
	
// }

/*instead of subtraction use division*/
// void solve(){
//   int A,B,N,b;
//   cin>>A>>B>>N;
//   vector<int> x(N);
//   for(int&a : x) cin >> a;
//   for(int&a : x){
//     cin>>b;
//     B -= (b+A-1)/A * a;
//   }
  
//   // hero can be dead after killing the last monster
//   // so B+m>0 where m is max health of all monsters
//   int m = 0;
//   for(int&a : x) m=max(m,a);
//   cout<<((B+m>0)?"YES\n":"NO\n");
// }

/*application of substr*/
// void solve(){
// int n;
// cin>>n;
// string s = "989";
// if(n<=3){cout<<s.substr(0,n)<<endl;return;}
//   cout<<s;
//   for(int i=3;i<n;i++)
//     cout<<(i-3)%10;
//   cout<<endl;
// }
/*value change concept, change += finalState-initialState*/
// bool check(int sum, int n) {
// 	// integer representation of sum / n >= 4.5
// 	return sum * 10 >= n * 45;
// }
// void solve(){
//   int n;
//   cin>>n;
//   std::vector<int> s(n);
//   int sum = 0;
//   for(int i=1;i<=n;i++)
//     cin>>s[i-1],sum+=s[i-1];
//   sort(s.begin(),s.end());
//   int cnt = 0;
//   while(!check(sum,n)){
//     int final = 5;
//     int orig = s[cnt];
//     sum += final-orig;
//     cnt++;
//   }
//   cout<<cnt<<endl;

// }
  



/** prev = next concept**/ 
// void solve(){
//   int n,k;
//   cin>>n>>k;
//   string s;cin>>s;
//   sort(s.begin(),s.end());
//   char last = 'a'-2;
//   int wt =0,len=0;
//   for(int i=0;i<n;i++){
//     if(s[i]>=last+2){
//       last = s[i];
//       wt +=s[i]-'a'+1;
//       len++;
//       if(len>=k){cout<<wt<<endl;return;}
//     }
//   }
//   cout<<"-1\n";
// }

  
/*function to check unique letters in string*/
// int check(string s,int n){
//   set<char> st;
//   for(int i=0;i<n;i++)st.insert(s[i]);
//   return st.size()==1;
// }

/*function to check palindromic string*/
// int isPalin(string s){
//   for(int i=0;i<s.size();i++)
//     if(s[i]!=s[s.size()-i-1])return false;
//   return true;
// }

/** reverse sort in vector and prefix sum(0-based) **/
// void solve(){
//   int n;
//   cin>>n;
//   vector<int> a(n);
//   for(auto &ele:a){
//     cin>>ele;
//   }
//   sort(a.rbegin(),a.rend());
//   int q;
//   cin>>q;
//   vector<int> pre(n+1,0);
//   pre[0]=0;
//   for(int i=1;i<=n;i++)
//     pre[i] = pre[i-1]+a[i-1];
//   while(q--){
//     int x;cin>>x;
//     int ans = pre[x-1];
//     int ans2 = pre[n]-pre[x];
  
//     cout<<ans+ans2<<endl;
//   }
  
// }

/*function to do run length encoding on a vector*/ 
// vector<pair<int,int>> rle(vector<int> nums){
//   vector<pair<int,int>> res;
//   if(nums.empty())return res;
//   int cnt = 1;
//   int curr = nums[0];
//   for(int i=1;i<nums.size();i++){
//     if(nums[i]==curr)cnt++;
//     else{
//       res.push_back({curr,cnt});
//       curr = nums[i];
//       cnt = 1;
//     }
//   }
//   res.push_back({curr,cnt});
//   return res;
// }

/** take out even and odd vector & use of accumulate to get the total sum **/
// void solve(){
//   int n;
//   cin>>n;
//   int sum = 0;
//   vector<int> even,odd;
//   for(int i=0;i<n;i++){
//     int x;cin>>x;
//     sum+=x;
//     if(x&1)odd.push_back(x);
//     else even.push_back(x);
//   }
//   sort(odd.rbegin(),odd.rend());
//   sort(even.rbegin(),even.rend());
  
//   int k = min(odd.size(),even.size());
//   int rv = sum;
//   rv -= accumulate(odd.begin(),odd.begin()+k,0);
//   rv -= accumulate(even.begin(),even.begin()+k,0);
//   if((odd.size())>k)rv-=odd[k];
//   if((even.size())>k)rv-=even[k];
//   cout<<rv<<endl;
// }

/** check filling of cross shaped tile entire across board **/ 
// void solve(){
//   int n; cin >> n;
 
//   char a[n][n];
//   for(int i=0; i<n; i++){
//       for(int j=0; j<n; j++){
//         cin >> a[i][j];
//       }
//   }
   
//   for(int i = 1; i < n - 1; i++) {
//       for(int j = 1; j < n - 1; j++) {
//         if(a[i][j] == '.' && a[i-1][j] == '.' && a[i+1][j] == '.' && a[i][j-1] == '.' && a[i][j+1] == '.') {
//             a[i][j] = a[i-1][j] = a[i+1][j] = a[i][j-1] = a[i][j+1] = '#';
//         }
//       }
//   }
   
//   for(int i = 0; i < n; i++) {
//       for(int j = 0; j < n; j++) {
//         if(a[i][j] == '.') {
//             cout << "NO\n";
//             return;
//         }
//       }
//   }
//   cout << "YES\n";
// }

/*direct printing using ternary operation*/
// void solve(){
//   int n,m;
//   cin>>n>>m;
//   cout<<(m?min(m,n-m):1)<<endl;
// }

/*find if any point is visited twice in a circle*/
// void solve(){
//   int n;cin>>n;
//   int p[n+1];
//   for(int i=1;i<=n;i++)cin>>p[i];
//   for(int i=1;i<=n;i++){
//     vector<int> seen(n+1,0);
//     // for(int i=1;i<=n;i++)cout<<seen[i]<<" ";
//     // cout<<endl;
    
//     int x  = p[i];
    
//     while(x!=i and !seen[x]){
      
//       seen[x]=1;
//       x = p[x];
//     }
//     cout<<x<<" ";
//   }
// }

/**function to check if splitting pile int x and 2x such that n=x/3+2x/3 reaches till m**/
// bool ok(int n,int m){
//   if(n==m)return 1;
//   else if(n%3!=0)return 0;
//   else {return(ok(n/3,m) or ok(2 * n / 3,m));}
// }


/*previous = next concept with swapping in a 2d grid*/
// void solve(){
//   int n,m;
//   cin>>n>>m;
//   char g[n + 7][m + 7];
// 	for (int i = 0; i < n; i++) {
// 		for (int j = 0; j < m; j++) {
// 			cin >> g[i][j];
// 		}
// 	}
// 	for (int j = 0; j < m; j++){
// 	  int last =  n-1;
// 	  for(int i=n-1;i>=0;i--){
// 	    if(g[i][j]=='o')last = i-1;
// 	    else if(g[i][j]=='*'){swap(g[i][j],g[last][j]);last--;}
// 	  }
// 	}
	
// 	for (int i = 0; i < n; i++) {
// 		for (int j = 0; j < m; j++) {
// 			cout<<g[i][j];
// 		}
// 		cout<<endl;
// 	}
// 	cout<<endl;
// }

/*DFS function to calculate number of leaf nodes(cnt[v]) for a particular vertice v*/
// vector<vector<int>> g;
// vector<int> cnt;
 
// void dfs(int v, int p) {
//     if (g[v].size() == 1 && g[v][0] == p) {
//         cnt[v] = 1;
//     } else {
//         for (auto u : g[v]) {
//             if (u != p) {
//                 dfs(u, v);
//                 cnt[v] += cnt[u];
//             }
//         }
//     }
// }

/*application of toggling state using set and keeping values between only 2 indices*/
/*also finding the clockwise and counterclockwise position from a starting position in a circle of length n*/
// void solve(){
//   int n, m, a; cin >> n >> m >> a;
//   set <int> q[2];
//   int ix = 0;
//   q[ix].insert(a);
//   while (m--) {
//     int x; char ch; cin>>x>>ch;
//     while(!q[ix].empty()){
//       int u = *(q[ix].begin());
//       q[ix].erase(u);
//       if(ch=='?' or ch=='0')
//         q[ix^1].insert((u+x-1)%n + 1);
      
//       if(ch=='?' or ch=='1')
//         q[ix^1].insert((u-x-1+n)%n + 1);
//     }
//     ix^=1;
//   }
//   cout<<q[ix].size()<<endl;
//   for(auto &x:q[ix])cout<<x<<" ";cout<<endl;
// }

/*function to find the size of cycle efficiently*/
// void solve(){
//   int n;
//   cin>>n;
//   int p[n+1];
//   for(int i=1;i<=n;i++)cin>>p[i];
  
//   vector<int> seen(n+1,0);
//   vector<int> ans(n+1,0);
  
//   for(int i=1;i<=n;i++){
//     if(seen[i])continue;
//     vector<int> cur;
    
    
//     while(!seen[i]){
//       cur.push_back(i);
//       seen[i]=1;
//       i = p[i];
//     }
//     for(auto el:cur)ans[el] = cur.size();
//   }
//   for(int i=1;i<=n;i++)cout<<ans[i]<<" ";
//   cout<<endl;
// }

/*DFS on 3d parallelopiped*/
// int n,m,k,ans,si,sj;
// char ch[15][15][15];
// void dfs(int x,int y,int z)
// {
// 	if (ch[x][y][z] == '#' || x>k || x<=0 || y>n || y<=0 || z>m || z<=0) return;
// 	ans++;
// 	ch[x][y][z]='#';
// 	dfs(x+1,y,z);
// 	dfs(x-1,y,z);
// 	dfs(x,y+1,z);
// 	dfs(x,y-1,z);
// 	dfs(x,y,z+1);
// 	dfs(x,y,z-1);
// }
// void solve(){
//   cin>>k>>n>>m;
// 	for(int i=1;i<=k;i++)
// 		for(int j=1;j<=n;j++)
// 			for(int k=1;k<=m;k++)
// 				cin>>ch[i][j][k];
// 	cin>>si>>sj;
// 	dfs(1,si,sj);
// 	cout<<ans<<endl;
// }


/*DFS
  function to find extra vals needed to be added in each subtree
  such that the overall value of all subtrees are equal starting from
  any vertice in a tree at the same level
*/
// int n;
// vector<int> a;

// int ans = 0;
// int dfs(int v){
//   int left = 2*v, right = 2*v + 1;
//   if(left>=(1<<(n+1)))return 0;
//   int L = a[left]+dfs(left);
//   int R = a[right]+dfs(right);
//   ans += abs(L-R);
//   return max(L,R);// returned from grandchild of left and right
// }
// void solve(){
//   cin>>n;
//   int m = (1<<(n+1));
//   a.assign(m,0);
//   for(int i=2;i<m;i++)cin>>a[i];
//   dfs(1);
//   cout<<ans<<endl;
// }

/* DFS function to find the maximum cost from every path in dfs using adjacency matrix */
// int ans;
// int g[110][110];
// int n;
  
// void dfs(int par,int node,int cost){
//   ans = max(ans,cost);
//   for(int nbr=0;nbr<n;nbr++){
//     if(nbr!=par and g[nbr][node])
//       dfs(node,nbr,cost+g[nbr][node]);
//   }
// }
// void solve(){
//   cin>>n;
//   for(int i=1;i<n;i++){
//     int u,v,c;cin>>u>>v>>c;
//     g[u][v] = c;g[v][u] = c;
//   }
//   dfs(-1,0,0);
//   cout<<ans<<endl;
// }
