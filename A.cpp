#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <set>
#include <utility>
#include <cstdio>
using namespace std;

class Transcoder {
public:
    virtual vector<vector<vector<int>>> init(int M, double e) = 0;
    virtual int decode(vector<vector<int>> H) = 0;
};

class TranscoderFixed: public Transcoder
{
public:
    vector<vector<vector<int>>> init(int M, double e)
    {
        int N = 4;
        vector<vector<vector<int>>> G(100, vector<vector<int>>(N, vector<int>(N)));
        return G;
    }

    int decode(vector<vector<int>> H)
    {
        return 0;
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
        }

        long long score = (long long)(1e9*pow(0.9, E)/N+.5);

        printf("M=%3d, e=%.2f, N=%3d, E=%3d, score=%8lld\n", M, e, N, E, score);

        total += score;
    }
    printf("total = %lld\n", total);
}

int main()
{
    TranscoderFixed fixed;
    Transcoder *transcoder = &fixed;

#ifdef WIN32
    test(transcoder);
#else
    int M;
    double e;
    cin>>M>>e;

    vector<vector<vector<int>>> G = transcoder->init(M, e);
    int N = G[0].size();
    cout<<N<<endl;
    for (int i=0; i<100; i++)
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
