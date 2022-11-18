#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <set>
#include <utility>
#include <cstdio>
#include <algorithm>
using namespace std;

class Transcoder {
public:
    virtual vector<vector<vector<int>>> init(int M, double e) = 0;
    virtual int decode(vector<vector<int>> H) = 0;
};

// 固定で0
class TranscoderFixed: public Transcoder
{
public:
    vector<vector<vector<int>>> init(int M, double e)
    {
        int N = 4;
        vector<vector<vector<int>>> G(M, vector<vector<int>>(N, vector<int>(N)));
        return G;
    }

    int decode(vector<vector<int>> H)
    {
        return 0;
    }
};

// 辺の本数にエンコード（N=100の辺を全て使う）
class TranscoderEdge: public Transcoder
{
    int M = 0;
    double e = 0.;
    int N = 0;
    int w = 0;
    vector<double> F;

    double comb(int a, int b)
    {
        return F[a]-F[b]-F[a-b];
    }

public:
    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;
        this->e = e;
        N = 100;

        w = (N*(N-1)/2-1)/(M-1);

        vector<vector<vector<int>>> G(M, vector<vector<int>>(N, vector<int>(N)));
        for (int i=0; i<M; i++)
        {
            int n = i*w;
            int c = 0;
            for (int j=0; j<N && c<n; j++)
                for (int k=0; k<j && c<n; k++)
                {
                    G[i][j][k] = G[i][k][j] = 1;
                    c++;
                }
        }

        F = vector<double>(N*N/2);
        F[0] = 0.;
        for (int i=1; i<N*N/2; i++)
            F[i] = log(i)+F[i-1];

        return G;
    }

    int decode(vector<vector<int>> H)
    {
        int c = 0;
        for (int i=0; i<N; i++)
            for (int j=0; j<i; j++)
                c += H[i][j];
        double e = max(0.01, this->e);

        double maxp = -1.;
        int maxi = 0;
        for (int i=0; i<M; i++)
        {
            int c0 = i*w;
            // c0 -> c となる確率
            double p = 0.;
            for (int j=0; j<=c0; j++)
            {
                int j2 = c-(c0-j);
                if (0<=j2 && j2<=N*(N-1)/2-c0)
                    p += exp(
                        comb(c0, j)+log(e)*j+log(1.-e)*(c0-j)+
                        comb(N*(N-1)/2-c0, j2)+log(e)*j2+log(1.-e)*(N*(N-1)/2-c0-j2));
            }

            if (p>maxp)
            {
                maxp = p;
                maxi = i;
            }
        }

        return maxi;
    }
};

vector<vector<string>> uniqueGraphs = {
    {}, {}, {}, {},
    {
        "000000", "000001", "000011", "000111", "001011", "001100", "001101", "001111",
        "011110", "011111", "111111",
    }, {
        "0000000000", "0000000001", "0000000011", "0000000111", "0000001011", "0000001100", "0000001101", "0000001111",
        "0000011110", "0000011111", "0000111111", "0001001011", "0001001100", "0001001101", "0001001111", "0001010110",
        "0001010111", "0001011110", "0001011111", "0001110100", "0001110101", "0001110111", "0001111111", "0011011110",
        "0011011111", "0011101011", "0011101100", "0011101101", "0011101111", "0011111110", "0011111111", "0111111011",
        "0111111111", "1111111111",
    }, {
        "000000000000000", "000000000000001", "000000000000011", "000000000000111", "000000000001011", "000000000001100", "000000000001101", "000000000001111",
        "000000000011110", "000000000011111", "000000000111111", "000000001001011", "000000001001100", "000000001001101", "000000001001111", "000000001010110",
        "000000001010111", "000000001011110", "000000001011111", "000000001110100", "000000001110101", "000000001110111", "000000001111111", "000000011011110",
        "000000011011111", "000000011101011", "000000011101100", "000000011101101", "000000011101111", "000000011111110", "000000011111111", "000000111111011",
        "000000111111111", "000001111111111", "000010001001011", "000010001001100", "000010001001101", "000010001001111", "000010001010100", "000010001010101",
        "000010001010110", "000010001010111", "000010001011110", "000010001011111", "000010001110100", "000010001110101", "000010001110111", "000010001111111",
        "000010010011110", "000010010011111", "000010010100000", "000010010100001", "000010010100011", "000010010100111", "000010010101010", "000010010101011",
        "000010010101100", "000010010101101", "000010010101110", "000010010101111", "000010010111110", "000010010111111", "000010011011110", "000010011011111",
        "000010011101011", "000010011101100", "000010011101101", "000010011101111", "000010011110100", "000010011110101", "000010011110110", "000010011110111",
        "000010011111110", "000010011111111", "000010110110000", "000010110110001", "000010110110011", "000010110110100", "000010110110101", "000010110110111",
        "000010110111001", "000010110111011", "000010110111100", "000010110111101", "000010110111111", "000010111111011", "000010111111100", "000010111111101",
        "000010111111111", "000011110110100", "000011110110101", "000011110110111", "000011110111111", "000011111111111", "000110011011110", "000110011011111",
        "000110011101010", "000110011101011", "000110011101100", "000110011101101", "000110011101110", "000110011101111", "000110011111110", "000110011111111",
        "000110101110000", "000110101110001", "000110101110011", "000110101110111", "000110101111000", "000110101111001", "000110101111011", "000110101111100",
        "000110101111101", "000110101111111", "000110111111000", "000110111111001", "000110111111010", "000110111111011", "000110111111110", "000110111111111",
        "000111100100001", "000111100100011", "000111100100111", "000111100101101", "000111100101111", "000111100111111", "000111101101011", "000111101101100",
        "000111101101101", "000111101101110", "000111101101111", "000111101110110", "000111101110111", "000111101111110", "000111101111111", "000111111111110",
        "000111111111111", "001110111111000", "001110111111001", "001110111111011", "001110111111111", "001111011011110", "001111011011111", "001111011101011",
        "001111011101100", "001111011101101", "001111011101111", "001111011111110", "001111011111111", "001111111111001", "001111111111011", "001111111111111",
        "011111111011110", "011111111011111", "011111111111111", "111111111111111",
    }
};

vector<vector<int>> normalize(vector<vector<int>> G)
{
    int N = (int)G.size();

    vector<vector<int>> R = G;

    vector<int> P;
    for (int i=0; i<N; i++)
        P.push_back(i);
    do
    {
        vector<vector<int>> G2(N, vector<int>(N));
        for (int i=0; i<N; i++)
            for (int j=0; j<i; j++)
                if (G[P[i]][P[j]]!=0)
                {
                    G2[i][j] = 1;
                    G2[j][i] = 1;
                }
        if (G2<R)
            R = G2;
    }
    while (next_permutation(P.begin(), P.end()));

    return R;
}

string to_string(vector<vector<int>> G)
{
    int N = (int)G.size();
    string g;
    for (int i=0; i<N; i++)
        for (int j=i+1; j<N; j++)
            g += "01"[G[i][j]];
    return g;
}

vector<vector<int>> from_string(string g)
{
    int N = 0;
    while (N*(N-1)/2<(int)g.size())
        N++;

    vector<vector<int>> G(N, vector<int>(N));
    int p = 0;
    for (int i=0; i<N; i++)
        for (int j=i+1; j<N; j++)
            G[i][j] = G[j][i] = g[p++]-'0';
    return G;
}

// 3個のクラスタのサイズにエンコード
class TranscoderCluster3: public Transcoder
{
    int M = 0;
    double e = 0.;
    vector<vector<int>> CS;

public:
    // パラメタを外部から与える
    bool tuning = false;
    // エラーが起こらないことを期待するモード
    bool exact = false;
    int N = 0;
    // クラスタサイズの間隔
    int CSint = 0;

    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;
        this->e = e;

        if (!tuning)
        {
            exact =
                e==0.0 ||
                abs(e-0.01)<1e-8 && M<=34 ||
                abs(e-0.02)<1e-8 && M<=11;

            int NT[11][11] = {
                {11, 16, 19, 22, 25, 27, 29, 31, 33, 35},
                {13, 16, 19, 22, 25, 27, 29, 32, 33, 35},
                {15, 18, 20, 23, 27, 27, 30, 32, 33, 36},
                {25, 23, 26, 27, 29, 31, 34, 35, 37, 35},
                {30, 32, 35, 32, 40, 36, 40, 39, 38, 45},
                {43, 41, 47, 47, 48, 44, 52, 53, 53, 55},
                {56, 59, 66, 63, 73, 64, 75, 66, 70, 67},
                {92, 89, 98, 96, 99, 92, 99, 99, 98, 97},
                {98, 100, 99, 96, 100, 98, 100, 100, 99, 99},
                {11, 99, 99, 96, 99, 98, 97, 96, 100, 100},
                {12, 16, 22, 23, 25, 27, 30, 31, 33, 37},
            };
            int ei = int(e*100+.5);
            N = NT[(ei+3)/4][M/10];

            CSint = 1;

            if (ei>=30 && M<=50)
            {
                N = 100;
                CSint = 3;
            }
        }

        if (exact)
        {
            N = 0;
            while ((int)uniqueGraphs[N].size()<M)
                N++;

            vector<vector<vector<int>>> G;
            for (int i=0; i<M; i++)
                G.push_back(from_string(uniqueGraphs[N][i]));
            return G;
        }

        // TODO: M種類のグラフを作れなくても、Nを小さくしたほうが良い場合もあるかも
        CS.clear();
        while (true)
        {
            for (int i=1; i<N; i++)
                for (int j=1; i+j<N; j++)
                    if (i<=j && j<=N-i-j && i%CSint==0 && j%CSint==0)
                        CS.push_back({i, j, N-i-j});
            // 最小のクラスタが大きいほど良い
            sort(CS.begin(), CS.end(), [](vector<int> &a, vector<int>&b) {
                return a[0]>b[0] || a[0]==b[0] && a[1]>b[1];
            });
            if (tuning && (int)CS.size()<M)
                return {};
            if ((int)CS.size()>=M)
                break;
            N++;
        }
        CS.resize(M);

        vector<vector<vector<int>>> G;
        for (int m=0; m<M; m++)
        {
            vector<vector<int>> g(N, vector<int>(N));
            for (int i=0; i<N; i++)
                for (int j=0; j<i; j++)
                    if (i<CS[m][0] && j<CS[m][0] ||
                        CS[m][0]<=i && i<CS[m][0]+CS[m][1] && CS[m][0]<=j && j<CS[m][0]+CS[m][1] ||
                        CS[m][0]+CS[m][1]<=i && CS[m][0]+CS[m][1]<=j)
                        g[i][j] = g[j][i] = 1;
            G.push_back(g);
        }

        return G;
    }

    int decode(vector<vector<int>> H)
    {
        if (exact)
        {
            string g = to_string(normalize(H));
            for (int i=0; i<M; i++)
                if (uniqueGraphs[N][i]==g)
                    return i;
            return 0;
        }

        vector<int> C(N);

        for (int iter=0; iter<16; iter++)
        {
            for (int i=0; i<N; i++)
            {
                int ms = -1;
                int mc = 0;
                for (int c=0; c<3; c++)
                {
                    int s = 0;
                    for (int j=0; j<N; j++)
                        if (j!=i)
                            if (H[i][j]!=0 && c==C[j] ||
                                H[i][j]==0 && c!=C[j])
                                s++;
                    if (s>ms)
                    {
                        ms = s;
                        mc = c;
                    }
                }
                C[i] = mc;
            }
        }

        vector<int> num(3);
        for (int i=0; i<N; i++)
            num[C[i]]++;
        sort(num.begin(), num.end());

        int ans = 0;
        for (int i=0; i<M; i++)
            if (abs(num[0]-CS[i][0])+abs(num[1]-CS[i][1]) <
                abs(num[0]-CS[ans][0])+abs(num[1]-CS[ans][1]))
                ans = i;
        return ans;
    }
};

void param_search()
{
    mt19937 rand;
    rand.seed(1234);

    TranscoderCluster3 transcoder;
    transcoder.tuning = true;

    for (int M=10; M<=100; M+=10)
        for (int ei=0; ei<=40; ei+=4)
        {
            double e = ei*.01;

            long long bestScore = 0;
            int bestN = 0;
            int bestCSint = 0;

            for (int CSint=1; CSint<=3; CSint+=2)
            for (int N=4; N<=101; N++)
            {
                transcoder.exact = N==101;
                transcoder.N = N;
                transcoder.CSint = CSint;

                vector<vector<vector<int>>> G = transcoder.init(M, e);
                long long score = 0;
                int Esum = 0;
                if (!G.empty())
                {
                    for (int i=0; i<10; i++)
                    {
                        int N = (int)G[0].size();
                        int E = 0;

                        for (int j=0; j<100; j++)
                        {
                            int s = rand()%M;
                            vector<vector<int>> H = G[s];
                            for (int k=0; k<N; k++)
                                for (int l=0; l<k; l++)
                                    if (int(rand()%100)<ei)
                                    {
                                        H[k][l] ^= 1;
                                        H[l][k] ^= 1;
                                    }
                            int ans = transcoder.decode(H);
                            if (ans!=s)
                                E++;
                        }
                        score += (long long)(1e9*pow(0.9, E)/N+.5);
                        Esum += E;
                    }
                }
                //printf("M = %d e = %.2f N = %d Csint = %d score = %lld Esum = %d\n", M, e, N, CSint, score, Esum);
                if (score>bestScore)
                {
                    bestScore = score;
                    bestN = N;
                    bestCSint = CSint;
                }
            }
            printf("%d %.2f %d %d %lld\n", M, e, bestN, bestCSint, bestScore);
        }
}

void test(Transcoder *transcoder)
{
    const int T = 50;

    long long total = 0LL;

    set<pair<int, int>> used;

    for (int i=0; i<T; i++)
    {
        mt19937 rand;
        rand.seed(i);

        int M;
        int ei;
        while (true)
        {
            M = rand()%91+10;
            ei = rand()%41;
            if (used.count({M, ei})==0)
                break;
        }
        used.insert({M, ei});
        double e = ei*.01;

        vector<vector<vector<int>>> G = transcoder->init(M, e);
        int N = (int)G[0].size();
        int E = 0;

        for (int j=0; j<100; j++)
        {
            int s = rand()%M;
            vector<vector<int>> H1 = G[s];
            for (int k=0; k<N; k++)
                for (int l=0; l<k; l++)
                    if (int(rand()%100)<ei)
                    {
                        H1[k][l] ^= 1;
                        H1[l][k] ^= 1;
                    }

            vector<int> P(N);
            for (int i=0; i<N; i++)
                P[i] = i;
            for (int i=0; i<N; i++)
                swap(P[i], P[rand()%(N-i)+i]);

            vector<vector<int>> H2(N, vector<int>(N));
            for (int k=0; k<N; k++)
                for (int l=0; l<k; l++)
                    if (H1[P[k]][P[l]]!=0)
                    {
                        H2[k][l] = 1;
                        H2[l][k] = 1;
                    }

            int ans = transcoder->decode(H2);
            if (ans!=s)
                E++;

            //printf("s=%d, ans=%d\n", s, ans);
        }

        long long score = (long long)(1e9*pow(0.9, E)/N+.5);

        printf("M=%3d, e=%.2f, N=%3d, E=%3d, score=%9lld\n", M, e, N, E, score);

        total += score;
    }
    printf("total = %lld\n", total);
}

void make_graph()
{
    for (int N=4; N<=6; N++)
    {
        cout<<"N="<<N<<endl;
        set<vector<vector<int>>> S;
        for (int b=0; b<1<<(N*(N-1)/2); b++)
        {
            vector<vector<int>> G(N, vector<int>(N));
            int c = 0;
            for (int i=0; i<N; i++)
                for (int j=0; j<i; j++)
                    if (b>>c++&1)
                        G[i][j] = G[j][i] = 1;
            S.insert(normalize(G));
        }
        for (vector<vector<int>> G: S)
            cout<<to_string(G)<<endl;
    }
}

int main()
{
    //param_search();
    //return 0;

    //make_graph();
    //return 0;

    //TranscoderFixed fixed;
    //TranscoderEdge edge;
    TranscoderCluster3 cluster3;
    Transcoder *transcoder = &cluster3;

#ifdef WIN32
    test(transcoder);
#else
    int M;
    double e;
    cin>>M>>e;

    vector<vector<vector<int>>> G = transcoder->init(M, e);
    int N = G[0].size();
    cout<<N<<endl;
    for (vector<vector<int>> g: G)
        cout<<to_string(g)<<endl;

    for (int k=0; k<100; k++)
    {
        string g;
        cin>>g;
        vector<vector<int>> H = from_string(g);
        int ans = transcoder->decode(H);
        cout<<ans<<endl;
    }
#endif
}
