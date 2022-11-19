#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <set>
#include <utility>
#include <cstdio>
#include <algorithm>
using namespace std;

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

vector<vector<int>> shuffle_graph(vector<vector<int>> G, mt19937 *rand)
{
    int N = (int)G.size();
    vector<int> P(N);
    for (int i=0; i<N; i++)
        P[i] = i;
    for (int i=0; i<N; i++)
        swap(P[i], P[(*rand)()%(N-i)+i]);

    vector<vector<int>> G2(N, vector<int>(N));
    for (int k=0; k<N; k++)
        for (int l=0; l<k; l++)
            if (G[P[k]][P[l]]!=0)
            {
                G2[k][l] = 1;
                G2[l][k] = 1;
            }
    return G2;
}

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

// 3個のクラスタのサイズにエンコード
class TranscoderCluster3: public Transcoder
{
    int M = 0;
    double e = 0.;
    vector<vector<int>> CS;

public:
    // パラメタを外部から与える
    bool tuning = false;
    int N = 0;
    // クラスタサイズの間隔
    int CSint = 0;

    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;
        this->e = e;

        if (!tuning)
        {
            int NT[14][18] = {
                {13, 15, 17, 19, 20, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35},
                {13, 15, 17, 19, 20, 22, 24, 24, 26, 27, 28, 29, 31, 31, 32, 33, 34, 35},
                {15, 17, 18, 19, 21, 22, 26, 26, 26, 27, 28, 31, 31, 32, 33, 34, 34, 35},
                {20, 20, 20, 20, 24, 24, 24, 26, 26, 28, 29, 32, 32, 32, 33, 36, 36, 36},
                {24, 27, 24, 27, 27, 28, 28, 32, 33, 34, 32, 34, 32, 35, 39, 38, 38, 39},
                {28, 30, 33, 32, 36, 35, 36, 38, 36, 33, 35, 37, 35, 39, 39, 41, 41, 44},
                {40, 40, 41, 44, 43, 42, 43, 41, 45, 51, 49, 42, 49, 51, 48, 44, 47, 53},
                {39, 53, 51, 55, 55, 56, 49, 50, 55, 61, 52, 63, 61, 61, 65, 62, 63, 57},
                {56, 56, 57, 69, 74, 78, 71, 71, 73, 73, 66, 71, 74, 74, 76, 76, 72, 75},
                {71, 75, 72, 73, 86, 90, 90, 94, 90, 93, 97, 100, 99, 99, 92, 99, 98, 98},
                {91, 92, 96, 91, 96, 98, 99, 99, 99, 100, 96, 100, 98, 99, 100, 100, 100, 100},
                {100, 100, 100, 100, 100, 100, 99, 99, 99, 99, 99, 99, 100, 100, 100, 99, 98, 98},
                {98, 100, 100, 99, 100, 99, 98, 98, 97, 97, 99, 99, 96, 100, 100, 100, 99, 98},
                {91, 93, 95, 97, 100, 24, 26, 24, 32, 27, 28, 30, 32, 33, 32, 35, 36, 37},
            };
            int CSintT[14][18] = {
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1},
                {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1},
                {5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1},
                {5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1},
                {5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            };
            int ei = int(e*100+.5);
            N = NT[ei/3][(M-10)/5];
            CSint = CSintT[ei/3][(M-10)/5];
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

// エラーが起こらないと仮定する
class TranscoderExact: public Transcoder
{
    int M = 0;
    int N = 0;
public:

    vector<vector<vector<int>>> init(int M, double e)
    {
        this->M = M;
        N = 0;
        while ((int)uniqueGraphs[N].size()<M)
            N++;

        vector<vector<vector<int>>> G;
        for (int i=0; i<M; i++)
            G.push_back(from_string(uniqueGraphs[N][i]));
        return G;
    }

    int decode(vector<vector<int>> H)
    {
        string g = to_string(normalize(H));
        for (int i=0; i<M; i++)
            if (uniqueGraphs[N][i]==g)
                return i;
        return 0;
    }
};

// トランスコーダーを使い分ける
class TranscoderIntegrate: public Transcoder
{
    TranscoderFixed fixed;
    TranscoderEdge edge;
    TranscoderCluster3 cluster3;
    TranscoderExact exact;

public:
    bool tuning = false;
    int id = 0;
    Transcoder *trans = nullptr;

    vector<vector<vector<int>>> init(int M, double e)
    {
        if (tuning)
        {
            switch (id) {
            case 0: trans = &fixed; break;
            case 1: trans = &edge; break;
            case 2: trans = &cluster3; break;
            case 3: trans = &exact; break;
            }
        }
        else
        {
            const int ID[41][91] = {
                {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
                {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
                {2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            };
            int id = ID[int(e*100+.5)][M-10];
            switch (id) {
            case 0: trans = &fixed; break;
            case 1: trans = &edge; break;
            case 2: trans = &cluster3; break;
            case 3: trans = &exact; break;
            }
        }

        return trans->init(M, e);
    }

    int decode(vector<vector<int>> H)
    {
        return trans->decode(H);
    }
};

void param_search()
{
    mt19937 rand;
    rand.seed(1234);

    auto test = [&](Transcoder *transcoder, int M, int ei, int n) -> long long
    {
        double e = ei*.01;
        vector<vector<vector<int>>> G = transcoder->init(M, e);
        long long score = 0;
        if (!G.empty())
        {
            for (int i=0; i<n; i++)
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
                    H = shuffle_graph(H, &rand);
                    int ans = transcoder->decode(H);
                    if (ans!=s)
                        E++;
                }
                score += (long long)(1e9*pow(0.9, E)/N+.5);
            }
        }
        return score;
    };

#if 0
    TranscoderCluster3 cluster3;
    cluster3.tuning = true;

    for (int ei=1; ei<=40; ei+=3)
    {
        vector<int> CS;

        for (int M=13; M<=100; M+=5)
        {
            long long bestScore = 0;
            int bestN = 0;
            int bestCSint = 0;

            for (int CSint=1; CSint<=5; CSint+=2)
            {
                long long preBestScore = 0;
                int preBestN = 0;
                for (int N=10; N<=100; N+=10)
                {
                    cluster3.N = N;
                    cluster3.CSint = CSint;

                    long long score = test(&cluster3, M, ei, 8);
                    if (score>preBestScore)
                    {
                        preBestScore = score;
                        preBestN = N;
                    }
                }
                for (int N=max(10, preBestN-9); N<=min(100, preBestN+9); N++)
                {
                    cluster3.N = N;
                    cluster3.CSint = CSint;

                    long long score = test(&cluster3, M, ei, 8);
                    if (score>bestScore)
                    {
                        bestScore = score;
                        bestN = N;
                        bestCSint = CSint;
                    }
                }
            }
            cout<<" "<<bestN;
            CS.push_back(bestCSint);
        }
        cout<<endl;
        for (int c: CS)
            cout<<" "<<c;
        cout<<endl;
    }

    cout<<endl;
#endif

#if 1
    TranscoderIntegrate transcoder;
    transcoder.tuning = true;

    for (int ei=40; ei<=40; ei++)
    {
        for (int M=10; M<=100; M++)
        {
            double e = ei*.01;

            long long bestScore = 0;
            int bestID = 0;

            for (int id=0; id<4; id++)
            {
                if (id==1)
                    if (!(ei>=35 || ei>=25 && 20<=M))
                        continue;

                transcoder.id = id;

                long long score = test(&transcoder, M, ei, 8);
                if (score>bestScore)
                {
                    bestScore = score;
                    bestID = id;
                }
            }
            cout<<" "<<bestID;
        }
        cout<<endl;
    }
#endif
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
            vector<vector<int>> H = G[s];
            for (int k=0; k<N; k++)
                for (int l=0; l<k; l++)
                    if (int(rand()%100)<ei)
                    {
                        H[k][l] ^= 1;
                        H[l][k] ^= 1;
                    }
            H = shuffle_graph(H, &rand);

            int ans = transcoder->decode(H);
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

    TranscoderFixed fixed;
    TranscoderEdge edge;
    TranscoderCluster3 cluster3;
    TranscoderExact exact;
    TranscoderIntegrate integrate;
    Transcoder *transcoder = &integrate;

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
