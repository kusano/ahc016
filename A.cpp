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

// 辺の本数にエンコード
class TranscoderEdgeNum: public Transcoder
{
    int N = 0;
    int M = 0;

public:
    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;

        N = 0;
        while (N*(N-1)/2+1<M)
            N++;

        vector<vector<vector<int>>> G(M, vector<vector<int>>(N, vector<int>(N)));
        for (int i=0; i<M; i++)
        {
            int c = 0;
            for (int j=0; j<N && c<i; j++)
                for (int k=0; k<j && c<i; k++)
                {
                    G[i][j][k] = G[i][k][j] = 1;
                    c++;
                }
        }

        return G;
    }

    int decode(vector<vector<int>> H)
    {
        int c = 0;
        for (int i=0; i<N; i++)
            for (int j=0; j<i; j++)
                c += H[i][j];
        return min(M-1, c);
    }
};

// 辺の本数にエンコード（N=100の辺を全て使う）
class TranscoderEdgeNum100: public Transcoder
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
        "0000000000000000", "0000000000010010", "0000000100010110", "0000001101010110", "0001000100011110", "0001001001001000", "0001001001011010", "0001001101011110",
        "0011001111001100", "0011001111011110", "0111101111011110",
    }, {
        "0000000000000000000000000", "0000000000000000000100010", "0000000000000010000100110", "0000000000000110010100110", "0000000001000010000101110", "0000000001000100010001000", "0000000001000100010101010", "0000000001000110010101110",
        "0000000011000110110001100", "0000000011000110110101110", "0000000111010110110101110", "0000100001000010000111110", "0000100001000100010011000", "0000100001000100010111010", "0000100001000110010111110", "0000100010000110110010100",
        "0000100010000110110110110", "0000100011000110110011100", "0000100011000110110111110", "0000100110010100110010000", "0000100110010100110110010", "0000100110010110110110110", "0000100111010110110111110", "0001100011000111110011100",
        "0001100011000111110111110", "0001100101010011000111110", "0001100101010101010011000", "0001100101010101010111010", "0001100101010111010111110", "0001100111010111110011100", "0001100111010111110111110", "0011100111110011100111110",
        "0011100111110111110111110", "0111110111110111110111110",
    }, {
        "000000000000000000000000000000000000", "000000000000000000000000000001000010", "000000000000000000000001000001000110", "000000000000000000000011000101000110", "000000000000000001000001000001001110", "000000000000000001000010000100001000", "000000000000000001000010000101001010", "000000000000000001000011000101001110",
        "000000000000000011000011001100001100", "000000000000000011000011001101001110", "000000000000000111001011001101001110", "000000000001000001000001000001011110", "000000000001000001000010000100011000", "000000000001000001000010000101011010", "000000000001000001000011000101011110", "000000000001000010000011001100010100",
        "000000000001000010000011001101010110", "000000000001000011000011001100011100", "000000000001000011000011001101011110", "000000000001000110001010001100010000", "000000000001000110001010001101010010", "000000000001000110001011001101010110", "000000000001000111001011001101011110", "000000000011000011000011011100011100",
        "000000000011000011000011011101011110", "000000000011000101001001010001011110", "000000000011000101001010010100011000", "000000000011000101001010010101011010", "000000000011000101001011010101011110", "000000000011000111001011011100011100", "000000000011000111001011011101011110", "000000000111000111011001011001011110",
        "000000000111000111011011011101011110", "000000001111010111011011011101011110", "000001000001000001000001000001111110", "000001000001000001000010000100111000", "000001000001000001000010000101111010", "000001000001000001000011000101111110", "000001000001000010000010001100110000", "000001000001000010000010001101110010",
        "000001000001000010000011001100110100", "000001000001000010000011001101110110", "000001000001000011000011001100111100", "000001000001000011000011001101111110", "000001000001000110001010001100110000", "000001000001000110001010001101110010", "000001000001000110001011001101110110", "000001000001000111001011001101111110",
        "000001000010000011000011011100101100", "000001000010000011000011011101101110", "000001000010000100001000010000100000", "000001000010000100001000010001100010", "000001000010000100001001010001100110", "000001000010000100001011010101100110", "000001000010000101001001010000101100", "000001000010000101001001010001101110",
        "000001000010000101001010010100101000", "000001000010000101001010010101101010", "000001000010000101001011010100101100", "000001000010000101001011010101101110", "000001000010000111001011011100101100", "000001000010000111001011011101101110", "000001000011000011000011011100111100", "000001000011000011000011011101111110",
        "000001000011000101001001010001111110", "000001000011000101001010010100111000", "000001000011000101001010010101111010", "000001000011000101001011010101111110", "000001000011000110001010011100110000", "000001000011000110001010011101110010", "000001000011000110001011011100110100", "000001000011000110001011011101110110",
        "000001000011000111001011011100111100", "000001000011000111001011011101111110", "000001000110000110011000011000100000", "000001000110000110011000011001100010", "000001000110000110011001011001100110", "000001000110000110011010011100100000", "000001000110000110011010011101100010", "000001000110000110011011011101100110",
        "000001000110000111011000011001101010", "000001000110000111011001011001101110", "000001000110000111011010011100101000", "000001000110000111011010011101101010", "000001000110000111011011011101101110", "000001000111000111011001011001111110", "000001000111000111011010011100111000", "000001000111000111011010011101111010",
        "000001000111000111011011011101111110", "000001001110010110011010011100100000", "000001001110010110011010011101100010", "000001001110010110011011011101100110", "000001001110010111011011011101101110", "000001001111010111011011011101111110", "000011000011000011000011111100111100", "000011000011000011000011111101111110",
        "000011000011000101001001110000111100", "000011000011000101001001110001111110", "000011000011000101001010110100111000", "000011000011000101001010110101111010", "000011000011000101001011110100111100", "000011000011000101001011110101111110", "000011000011000111001011111100111100", "000011000011000111001011111101111110",
        "000011000101000110011000101000110000", "000011000101000110011000101001110010", "000011000101000110011001101001110110", "000011000101000110011011101101110110", "000011000101000111011000101000111000", "000011000101000111011000101001111010", "000011000101000111011001101001111110", "000011000101000111011010101100111000",
        "000011000101000111011010101101111010", "000011000101000111011011101101111110", "000011000111000111011000111000111000", "000011000111000111011000111001111010", "000011000111000111011001111000111100", "000011000111000111011001111001111110", "000011000111000111011011111100111100", "000011000111000111011011111101111110",
        "000011001100010100011000100001100010", "000011001100010100011001100001100110", "000011001100010100011011100101100110", "000011001100010101011010100101101010", "000011001100010101011011100101101110", "000011001100010111011011101101101110", "000011001101010101011001100001111110", "000011001101010101011010100100111000",
        "000011001101010101011010100101111010", "000011001101010101011011100100111100", "000011001101010101011011100101111110", "000011001101010110011011101100110100", "000011001101010110011011101101110110", "000011001101010111011011101100111100", "000011001101010111011011101101111110", "000011001111010111011011111100111100",
        "000011001111010111011011111101111110", "000111000111000111111000111000111000", "000111000111000111111000111001111010", "000111000111000111111001111001111110", "000111000111000111111011111101111110", "000111001011010011100011111100111100", "000111001011010011100011111101111110", "000111001011010101101001110001111110",
        "000111001011010101101010110100111000", "000111001011010101101010110101111010", "000111001011010101101011110101111110", "000111001011010111101011111100111100", "000111001011010111101011111101111110", "000111001111010111111000111001111010", "000111001111010111111001111001111110", "000111001111010111111011111101111110",
        "001111001111110011110011111100111100", "001111001111110011110011111101111110", "001111001111110111111011111101111110", "011111101111110111111011111101111110",
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

            vector<vector<vector<int>>> G(M, vector<vector<int>>(N, vector<int>(N)));
            for (int m=0; m<M; m++)
                for (int i=0; i<N; i++)
                    for (int j=0; j<N; j++)
                        G[m][i][j] = uniqueGraphs[N][m][i*N+j]-'0';
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
            H = normalize(H);
            string g;
            for (int i=0; i<N; i++)
                for (int j=0; j<N; j++)
                    g += "01"[H[i][j]];
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

// 2個のクラスタのサイズにエンコード
class TranscoderCluster2: public Transcoder
{
    int M = 0;
    double e = 0.;
    int N = 0;
    vector<vector<int>> CS;

public:
    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;
        this->e = e;
        N = 100;

        CS.clear();
        for (int i=1; i<N; i++)
            for (int j=1; i+j<=N; j++)
                if (i<=j)
                    CS.push_back({i, j});
        // 小さいクラスタが大きいほど良い
        sort(CS.begin(), CS.end(), [](vector<int> &a, vector<int>&b) {
            return a[0]>b[0] || a[0]==b[0] && a[1]>b[1];
        });
        CS.resize(M);

        vector<vector<vector<int>>> G;
        for (int m=0; m<M; m++)
        {
            vector<vector<int>> g(N, vector<int>(N));
            for (int i=0; i<N; i++)
                for (int j=0; j<i; j++)
                    if (i<CS[m][0] && j<CS[m][0] ||
                        CS[m][0]<=i && i<CS[m][0]+CS[m][1] && CS[m][0]<=j && j<CS[m][0]+CS[m][1])
                        g[i][j] = g[j][i] = 1;
            G.push_back(g);
        }

        return G;
    }

    int decode(vector<vector<int>> H)
    {
        vector<int> C(N, 2);

        for (int iter=0; iter<16; iter++)
        {
            for (int i=0; i<N; i++)
            {
                int ms = -1;
                int mc = 0;
                for (int c=0; c<2; c++)
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

        for (int iter=0; iter<16; iter++)
        {
            vector<int> num(3);
            for (int i=0; i<N; i++)
                num[C[i]]++;

            for (int i=0; i<N; i++)
            {
                int c[2] = {};
                for (int j=0; j<N; j++)
                    if (H[i][j]!=0 && C[j]<=1)
                        c[C[j]]++;
                if (c[0]>c[1] && c[0]*2>=num[0]-(C[i]==0?1:0))
                    C[i] = 0;
                else if (c[1]>c[0] && c[1]*2>=num[1]-(C[i]==1?1:0))
                    C[i] = 1;
                else
                    C[i] = 2;
            }
        }

        vector<int> num(2);
        for (int i=0; i<N; i++)
            if (C[i]<2)
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
        {
            string s;
            for (int i=0; i<N; i++)
                for (int j=0; j<N; j++)
                    s += "01"[G[i][j]];
            cout<<s<<endl;
        }
    }
}

int main()
{
    //param_search();
    //return 0;

    //make_graph();
    //return 0;

    //TranscoderFixed fixed;
    //TranscoderEdgeNum edge;
    //TranscoderEdgeNum100 edge100;
    TranscoderCluster3 cluster3;
    //TranscoderCluster2 cluster2;
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
    for (int i=0; i<M; i++)
    {
        string g;
        for (int j=0; j<N; j++)
            for (int k=j+1; k<N; k++)
                g += "01"[G[i][j][k]];
        cout<<g<<endl;
    }

    for (int k=0; k<100; k++)
    {
        string g;
        cin>>g;
        int p = 0;
        vector<vector<int>> H(N, vector<int>(N));
        for (int i=0; i<N; i++)
            for (int j=i+1; j<N; j++)
                H[i][j] = H[j][i] = g[p++]-'0';
        int ans = transcoder->decode(H);
        cout<<ans<<endl;
    }
#endif
}
