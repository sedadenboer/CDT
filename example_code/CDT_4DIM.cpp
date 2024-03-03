#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//Constants
#define DIM 4 //3
#define DIMEN DIM+1
#define MSUB 10 //6
#define MXRNDSM 10000000 //65536
#define MAXPNTS 10000000 //65536

typedef unsigned int vertex;

struct ssimplex {
    struct ssimplex *next[DIMEN];
    vertex point[DIMEN];
    vertex nxtpnt[DIMEN];
    unsigned int flag;
};
typedef struct ssimplex simplex;

static int numsubsim[DIMEN];
static int subsim[DIMEN][MSUB][DIMEN];
static int deltan[DIMEN][DIMEN];
int nsims[DIMEN];
int simalloc;

simplex *lattice;
vertex pntalloc;

static int combi(int n, int k) {
    if ((k < 0) || (k > n))
        return 0;
    if(n==k)
        return 1;

    int i, nf = 1, kf = 0, nmkf = 0;
    for(i = 1; i < n +1; i++) {
        nf *= i;
        if (i == k) kf = nf;
        if (i == (n-k)) nmkf = nf;
    }

    return nf / (kf*nmkf);
}

unsigned int randbb(void){
    return rand();
}

void initmove(void)
{
    int i;
    numsubsim[0] = 1;

    for(i=1; i<DIMEN; i++)
        numsubsim[i]=numsubsim[i-1]*(DIMEN+1-i) / i;

    memset(subsim, -1, sizeof(subsim));
    for(i=0; i<DIMEN; i++) {
        int j;
        int nowsubs [DIMEN];
        for (j=0; j<i; j++)
            nowsubs [j] = j;
        for (j=i; j<DIMEN; j++)
            nowsubs [j] = -1;
        for (j=0; j<numsubsim[i]; j++) {
            int k;
            int neighbor = i-1;
            int ssind = 0, nsind = i;
            memcpy(subsim[i][j], nowsubs, sizeof(nowsubs));
            for (k=0; k<DIMEN; k++)
                if (subsim[i][j][ssind] == k)
                    ssind++;
                else
                    subsim[i][j][nsind++] = k;
                while (++nowsubs[neighbor] == DIMEN + 1 + neighbor - i)
                    neighbor--;
                for (k=neighbor+1; k<i; k++)
                    nowsubs[k] = nowsubs[k-1] + 1;
        }
    }

    for (i=0; i<DIMEN; i++) {
        int j;
        for (j=0; j<DIMEN; j++)
            deltan[i][j] = combi(DIMEN-i, DIMEN-j) - combi(i+1, DIMEN-j);
    }
}

void makesmall(void) {
    int i;
    simalloc = DIMEN + 1;
    lattice = (simplex *)sbrk(simalloc * sizeof(simplex));
    pntalloc = DIMEN + 1;
    for (i=0; i<DIMEN+1; i++) {
        int j;
        int k = 0;
        for (j=0; j<DIMEN+1; j++) {
            if (j == i)
                continue;
            lattice[i].next[k] = lattice + j;
            lattice[i].point[k] = j;
            k++;
        }
     for (j=0; j<DIMEN; j++)
        lattice[i].nxtpnt[j] = i;
        lattice[i].flag = 0;
   }
   for (i=0; i<DIMEN; i++)
     nsims[i] = combi(DIMEN+1, i+1);
}

int ndmove(int movetype)
{
    int subs;
    simplex *base;
    simplex *moving [DIMEN+1];
    vertex mpoint [DIMEN+1];
    int i;
    int *heresub;
    vertex newpoint;
    static unsigned int recount = 0;
    static simplex *first = 0;
    static int frstpnt = 0;
    #define MXFRPNT 81920
    static vertex freepnts [MXFRPNT];

    do {
        base = lattice + (randbb() % simalloc);
    } while (!base->next[0]);
    subs = randbb() % numsubsim[movetype];
    heresub = subsim[movetype][subs];

    if (movetype) {
        int frst, last;
        static simplex *list[MXRNDSM];
        static char repnt[MAXPNTS];
        newpoint = base->nxtpnt[heresub[0]];

        for (i=1; i<movetype; i++)
         if (base->nxtpnt[heresub[i]] != newpoint)
            return 1;

        for (; i<DIMEN; i++)
          if (base->nxtpnt[heresub[i]] == newpoint)
            return 2;


        if (++recount == 0) {
            for (i=0; i<simalloc; i++)
                lattice[i].flag = 0;
                recount = 1;
        }
        base->flag = recount;

        frst = 0;
        last = 0;

        for (i = 0; i < movetype; i++)
            repnt[base->point[heresub[i]]] = 1;

        for (; i < DIMEN; i++)
            (list[last++] = base->next[heresub[i]])->flag = recount;

        do {
            simplex *thissim = list[frst++];
            for (i = 0; i < DIMEN; i++)
                if (thissim->point[i] == newpoint) {
                    for (i = 0; i < movetype; i++)
                        repnt[base->point[heresub[i]]] = 0;
                    return 3;
                }
            for (i = 0; i < DIMEN; i++) {
                if (repnt[thissim->point[i]] ||
                    thissim->next[i]->flag == recount)
                    continue;
                (list[last++] = thissim->next[i])->flag = recount;
            }
        } while (frst < last);

        for (i = 0; i < movetype; i++)
            repnt[base->point[heresub[i]]] = 0;
    } else {
        if (frstpnt)
            newpoint = freepnts[--frstpnt];
        else {
            newpoint = pntalloc++;
            if (!pntalloc) {
                fprintf(stderr,"Too many points for vertex data type!\n");
                exit(1);
            }
        }
    }
    moving[0] = base;
    for (i=1; i<movetype+1; i++)
        moving[i] = base->next[heresub[i-1]];
    for (; i<DIMEN+1; i++) {
        if (!first) {
            simalloc++;
            brk((char *)(lattice + simalloc));
            first = lattice + simalloc - 1;
            first->next[0] = 0;
            first->next[1] = 0;
            first->flag = 0;
        }
        moving[i] = first;
        first = first->next[1];
    }
    mpoint[0] = newpoint;
    for (i=1; i<DIMEN+1; i++)
        mpoint[i] = base->point[heresub[i-1]];
    for (i=movetype+1; i<DIMEN+1; i++) {
        int j;
        int k = 0;
        for (j=0; j<DIMEN+1; j++) {
            if (j == i)
                continue;
            if (j > movetype) {
                moving[i]->next[k] = moving[j];
                moving[i]->nxtpnt[k] = mpoint[i];
            } else {
                simplex *oldsim = moving[j];
                int oldtonon;
                simplex *nonsim;
                int nontoold;
                for (oldtonon=0; oldtonon<DIMEN; oldtonon++)
                    if (oldsim->point[oldtonon] == mpoint[i])
                        break;
                nonsim = oldsim->next[oldtonon];
                for (nontoold=0; nontoold<DIMEN; nontoold++)
                    if (nonsim->next[nontoold] == oldsim)
                        break;
                moving[i]->next[k] = nonsim;
                moving[i]->nxtpnt[k] = nonsim->point[nontoold];
                nonsim->next[nontoold] = moving[i];
                nonsim->nxtpnt[nontoold] = mpoint[j];
            }
            moving[i]->point[k] = mpoint[j];
            k++;
        }
    }
    for (i=0; i<movetype+1; i++) {
        moving[i]->next[0] = 0;
        moving[i]->next[1] = first;
        first = moving[i];
    }
    if (movetype == DIMEN-1) {
        if (frstpnt == MXFRPNT) {
            fprintf(stderr,"Size of freepnt too small!\n");
            exit(1);
        }
        freepnts[frstpnt++] = mpoint[DIMEN];
    }
    for (i=0; i<DIMEN; i++)
        nsims[i] += deltan[movetype][i];
    return 0;
}

void printInfo(void){
    printf("Simalloc: %d\n", simalloc);
    int i, j;
    for(i = 0; i< simalloc; i++){
        for(j=0; j < DIMEN; j++)
            printf("%d\t", (int)lattice[i].point[j]);
        printf("\n");
    }
}

/*
 *Main entry point to the program
 */
int main(int argc, char** argv) {
    int res, i, good, steps;
    int moves[DIMEN][2];

    if(argc < 2) {
        printf("Usage: %s [number of moves]\n", argv[0]);
        printf("Defaulting to 100\n");
        steps = 100;
    }
    else
        steps = atoi(argv[1]);

    printf("Taking %d steps...\n", steps);

    for(i=0; i<DIMEN; i++) {
        moves[i][0] = 0;
        moves[i][1] = 0;
    }

    srand( (unsigned int)time( NULL ) );

    initmove();
    makesmall();

    //printInfo();

    int toMove[3][2] = {0, 4,
                        1, 3,
                        2, 2};

    good = 0;
    for(i=0; i < steps; i++)
    {
        int move_group = rand() % 3;
        int move = rand() % 2;

        if(!toMove[move_group][move]) {
            //if(rand() % ((pntalloc / 10) + 1))
            if(nsims[0] > 32000)
                move = 1;
        }

        res = ndmove(toMove[move_group][move]);

        moves[toMove[move_group][move]][0]++;

        if(!res){
            good++;
            moves[toMove[move_group][move]][1]++;
        }
    }
    printf("\nNumber of good moves: %d - %f%\n", good, ((float)good/(float)steps)*100);

    for(i=0; i<DIMEN; i++)
        printf("Move %d = %d / %d\n", i, moves[i][1],moves[i][0]);

    for(i=0; i < DIMEN; i++)
        printf("N_%d = %d \n", i, nsims[i]);



    //if (res != 0)
    //    printf("Test %d failed!\n", res);

    //printInfo();

    return (EXIT_SUCCESS);
}
