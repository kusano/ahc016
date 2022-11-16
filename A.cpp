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

void test(Transcoder *transcoder)
{
    const int T = 50;

    mt19937 rand;
    rand.seed(1234);

    long long total = 0LL;

    set<pair<int, int>> used;

    for (int i=0; i<T; i++)
    {
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

        printf("M=%3d, e=%.2f, N=%3d, E=%3d, score=%8lld\n", M, e, N, E, score);

        total += score;
    }
    printf("total = %lld\n", total);
}

int main()
{
    //TranscoderFixed fixed;
    //TranscoderEdgeNum edge;
    TranscoderEdgeNum100 edge100;
    Transcoder *transcoder = &edge100;

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
