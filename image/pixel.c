
#include "pixel.h"

void pixel_avg(pixel_t *out, const pixel_t **in, const int *weight, int count, int mask) {
        wpixel_t p={};
        int totalweight=0;
        for(int i=0; i<count; i++) {
                wpixel_add(&p,in[i],weight[i],mask);
                totalweight += weight[i];
        }
        wpixel_div(&p,totalweight,mask);

        for(int c=0; c<3; c++) {
           if (!(mask&(1<<c))) {
             out->d[c] = p.wd[c];
           }
        }
}
void pixel_med(pixel_t *out, const pixel_t **in, const int *weight, int count, int mask) {
        wpixel_t p={};
        int totalweight=0;
        for(int i=0; i<count; i++) {
                wpixel_add(&p,in[i],weight[i],mask);
                totalweight += weight[i];
        }
        wpixel_div(&p,totalweight,mask);

        for(int c=0; c<3; c++) {
           if (!(mask&(1<<c))) {
             int min_dist=100000;
             int min_pos=-1;
             for(int i=0; i<count; i++) {
                int diff=abs(p.wd[c]-in[i]->d[c]);
                if (min_pos==-1 || diff<min_dist) {min_dist=diff; min_pos=i;}
             }
             out->d[c] = in[min_pos]->d[c];
           }
        }

}
