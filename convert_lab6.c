#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEIGHT        160
#define WIDTH         120
#define HEIGHT_START  30
#define WIDTH_START   65
#define HEIGHT_END    65
#define WIDTH_END     100
#define HEIGHT_START2 3
#define WIDTH_START2  35
#define HEIGHT_END2   50
#define WIDTH_END2    90


void ipr_load_ppm(unsigned char image[][WIDTH][3], const char path[]);
void ipr_save_ppm(unsigned char image[][WIDTH][3], const char path[]);
void gradation_plate_range(unsigned char rgb[][WIDTH][3], unsigned char rgb2[][WIDTH][3], int x, int y, int z, int w);
void rgb3lab(unsigned char rgb[][WIDTH][3], double lab[][WIDTH][3]);
void convert_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3], double lab3[][WIDTH][3]);
void lab3rgb(double lab[][WIDTH][3], unsigned char rgb[][WIDTH][3]);
void gradation_plate_lab(double lab[][WIDTH][3], double*, double*, double*, double*, int x, int y, int z, int w);
void gradation_plate_lab2(double lab[][WIDTH][3], double*, double*, double*, double*, int x, int y, int z, int w);
void gradation_plate_move_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3], double*, double*, double*, double*, double*, double*, double*, double*, double lab3[][WIDTH][3], double lab4[][WIDTH][3]);
void gradation_plate_move_lab2(double lab[][WIDTH][3], double lab2[][WIDTH][3], double*, double*, double*, double*, double*, double*, double*, double*, double lab3[][WIDTH][3], double lab4[][WIDTH][3]);

int main(int argc, char *argv[])
{
    unsigned char src_img[HEIGHT][WIDTH][3];
    unsigned char src_img2[HEIGHT][WIDTH][3];
    int dst_img[HEIGHT][WIDTH][3];
    int dst_img2[HEIGHT][WIDTH][3];    
    unsigned char img[HEIGHT][WIDTH][3];
    unsigned char img2[HEIGHT][WIDTH][3];
    int plate[HEIGHT][WIDTH][3];
    int plate2[HEIGHT][WIDTH][3];
    double hensa1;
    double hensa2;
    double R1;
    double R2;
    int m,n;
    unsigned char range[HEIGHT][WIDTH][3];
    unsigned char range2[HEIGHT][WIDTH][3];
    double a1,b1,c1,a2,b2,c2;
    double img_lab[HEIGHT][WIDTH][3];
    double img_lab2[HEIGHT][WIDTH][3];
    double img_lab3[HEIGHT][WIDTH][3];
    unsigned char rgb_lab[HEIGHT][WIDTH][3];
    FILE *fp;
    FILE *fp2;
    double a3,b3,c3,a4,b4,c4,a5,b5,c5,a6,b6,c6;
    double hensa3=0;
    double hensa4;
    double hensa5=0;
    double hensa6;
    double lab2[HEIGHT][WIDTH][3];
    unsigned char rgb_lab2[HEIGHT][WIDTH][3];
    FILE *fp3;
    double lab_plate[HEIGHT][WIDTH][3];
    double lab_plate2[HEIGHT][WIDTH][3];
    FILE *fp4;

    if (argc != 7) {
      fprintf(stderr, "Usage: %s source destination\n", argv[0]);
        exit(1);
    }

    ipr_load_ppm(src_img, argv[1]);
    ipr_load_ppm(src_img2, argv[2]);

    gradation_plate_range(src_img, range, HEIGHT_START, WIDTH_START, HEIGHT_END, WIDTH_END);
    gradation_plate_range(src_img2, range2,  HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);

    ipr_save_ppm(range, argv[3]);
    ipr_save_ppm(range2, argv[4]);

    rgb3lab(src_img, img_lab);
    rgb3lab(src_img2, img_lab2);

    fp=fopen("lab.txt","w");
    for(m=HEIGHT_START;m<HEIGHT_END;m++){
      for(n=WIDTH_START;n<WIDTH_END;n++){
	fprintf(fp,"%lf %lf %lf\n",img_lab[m][n][1],img_lab[m][n][2],img_lab[m][n][0]);
      }
    }
    fclose(fp);

    fp2=fopen("lab2.txt","w");
    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
      for(n=WIDTH_START2;n<WIDTH_END2;n++){
	fprintf(fp2,"%lf %lf %lf\n",img_lab2[m][n][1],img_lab2[m][n][2],img_lab2[m][n][0]);
      }
    }
    fclose(fp2);

    for(m=0;m<HEIGHT;m++){
        for(n=0;n<WIDTH;n++){
            lab_plate[m][n][0]=img_lab[m][n][0];
            lab_plate[m][n][1]=img_lab[m][n][1];
            lab_plate[m][n][2]=img_lab[m][n][2];
            lab_plate2[m][n][0]=img_lab2[m][n][0];
            lab_plate2[m][n][1]=img_lab2[m][n][1];
            lab_plate2[m][n][2]=img_lab2[m][n][2];
        }
    }

    convert_lab(img_lab, img_lab2, img_lab3);
    lab3rgb(img_lab3, rgb_lab);

    ipr_save_ppm(rgb_lab, argv[5]);

    gradation_plate_lab(lab_plate, &a3, &b3, &hensa3, &hensa4, HEIGHT_START, WIDTH_START, HEIGHT_END, WIDTH_END);
    if(hensa3==0){
      hensa6=1;
    gradation_plate_lab(lab_plate2, &a4, &b4, &hensa5, &hensa6,  HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);
    }else{
      gradation_plate_lab2(lab_plate2, &a4, &b4, &hensa5, &hensa6,  HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);
    }

    if(hensa3==0){
      gradation_plate_move_lab(lab_plate, lab2, &a3, &b3, &a4, &b4, &hensa3, &hensa5,  &hensa4, &hensa6, img_lab3, lab_plate2);
    }else{
      gradation_plate_move_lab2(lab_plate, lab2, &a3, &b3, &a4, &b4, &hensa3, &hensa5,  &hensa4, &hensa6, img_lab3, lab_plate2);
    }

    lab3rgb(lab2, rgb_lab2);

    fp3=fopen("lab3.txt","w");
    for(m=HEIGHT_START;m<HEIGHT_END;m++){
        for(n=WIDTH_START;n<WIDTH_END;n++){
            fprintf(fp3,"%lf %lf %lf\n",lab2[m][n][1],lab2[m][n][2],lab2[m][n][0]);
      }
    }
    fclose(fp3);

    fp4=fopen("lab4.txt","w");
    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
      for(n=WIDTH_START2;n<WIDTH_END2;n++){
	fprintf(fp4,"%lf %lf %lf\n",img_lab3[m][n][1],img_lab3[m][n][2],img_lab3[m][n][0]);
      }
    }
    fclose(fp4);

    ipr_save_ppm(rgb_lab2, argv[6]);
    return 0;
}


void ipr_load_ppm(unsigned char image[][WIDTH][3], const char path[])
{
    char magic_number[2];
    int width, height;
    int max_intensity;
    FILE *fp;
    
    fp = fopen(path, "rb");
    if (fp == NULL) {
        fprintf(stderr, "%s が開けませんでした．\n", path);
        exit(1);
    }
    
    fscanf(fp, "%c%c", &magic_number[0], &magic_number[1]);
    if (magic_number[0] != 'P' || magic_number[1] != '6') {
        fprintf(stderr, "%s はバイナリ型 PPM ではありません．\n", path);
        fclose(fp);
        exit(1);
    }
    
    fscanf(fp, "%d %d", &width, &height);
    if (width != WIDTH || height != HEIGHT) {
        fprintf(stderr, "画像のサイズが異なります．\n");
        fprintf(stderr, "想定サイズ：WIDTH = %d, HEIGHT = %d\n", WIDTH, HEIGHT);
        fprintf(stderr, "実サイズ：  width = %d, height = %d\n", width, height);
        fclose(fp);
        exit(1);
    }
    
    fscanf(fp, "%d", &max_intensity);
    if (max_intensity != 255) {
        fprintf(stderr, "最大階調値が不正な値です（%d）．\n", max_intensity);
        fclose(fp);
        exit(1);
    }
    
    fgetc(fp); 

    fread(image, sizeof(unsigned char), HEIGHT * WIDTH * 3, fp);
    
    fclose(fp);
}

void ipr_save_ppm(unsigned char image[][WIDTH][3], const char path[])
{
    FILE *fp;
    
    fp = fopen(path, "wb");
    if (fp == NULL) {
        fprintf(stderr, "%s が開けませんでした．\n", path);
        exit(1);
    }
    
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", WIDTH, HEIGHT);
    fprintf(fp, "255\n");
    fwrite(image, sizeof(unsigned char), HEIGHT * WIDTH * 3, fp);
    
    fclose(fp);
}

void gradation_plate_range(unsigned char rgb[][WIDTH][3], unsigned char rgb2[][WIDTH][3], int x, int y, int z, int w)
{
  int m,n;

  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      rgb2[m][n][0]=rgb[m][n][0];
      rgb2[m][n][1]=rgb[m][n][1];
      rgb2[m][n][2]=rgb[m][n][2];
    }
  }

  for(m=x;m<x+1;m++){
    for(n=y;n<w;n++){
      rgb2[m][n][0]=1;
      rgb2[m][n][1]=255;
      rgb2[m][n][2]=255;
    }
  }

  for(m=z;m<z+1;m++){
    for(n=y;n<w;n++){
      rgb2[m][n][0]=1;
      rgb2[m][n][1]=255;
      rgb2[m][n][2]=255;
    }
  }

  for(m=x;m<z;m++){
    for(n=y;n<y+1;n++){
      rgb2[m][n][0]=1;
      rgb2[m][n][1]=255;
      rgb2[m][n][2]=255;
    }
  }

  for(m=x;m<z;m++){
    for(n=w;n<w+1;n++){
      rgb2[m][n][0]=1;
      rgb2[m][n][1]=255;
      rgb2[m][n][2]=255;
    }
  }
}

void rgb3lab(unsigned char rgb[][WIDTH][3], double lab[][WIDTH][3])
{
  int m,n;
  double rgb2[HEIGHT][WIDTH][3];
  double xyz[HEIGHT][WIDTH][3];

  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      rgb2[m][n][0]=(double)rgb[m][n][0]/255;
      rgb2[m][n][1]=(double)rgb[m][n][1]/255;
      rgb2[m][n][2]=(double)rgb[m][n][2]/255;
      //printf("%lf %lf %lf\n",rgb2[m][n][0],rgb2[m][n][1],rgb2[m][n][2]);

      if(rgb2[m][n][0]<=0.04045){
	rgb2[m][n][0]/=12.92;
      }else{
	rgb2[m][n][0]=pow((rgb2[m][n][0]+0.055)/1.055,2.4);
      }
      if(rgb2[m][n][1]<=0.04045){
	rgb2[m][n][1]/=12.92;
      }else{
	rgb2[m][n][1]=pow((rgb2[m][n][1]+0.055)/1.055,2.4);
      }
      if(rgb2[m][n][2]<=0.04045){
	rgb2[m][n][2]/=12.92;
      }else{
	rgb2[m][n][2]=pow((rgb2[m][n][2]+0.055)/1.055,2.4);
      }

      xyz[m][n][0]=0.436041*rgb2[m][n][0]+0.385113*rgb2[m][n][1]+0.143046*rgb2[m][n][2];
      xyz[m][n][1]=0.222485*rgb2[m][n][0]+0.716905*rgb2[m][n][1]+0.060610*rgb2[m][n][2];
      xyz[m][n][2]=0.013920*rgb2[m][n][0]+0.097067*rgb2[m][n][1]+0.713913*rgb2[m][n][2];  
      //printf("%lf %lf %lf\n",xyz[m][n][0],xyz[m][n][1],xyz[m][n][2]);
      if(xyz[m][n][0]>0.008856){
	xyz[m][n][0]=pow(xyz[m][n][0],1.0/3.0);
      }else{
	xyz[m][n][0]=(pow(29.0/3.0,3)*xyz[m][n][0]+16)/116;
      }
      if(xyz[m][n][1]>0.008856){
	xyz[m][n][1]=pow(xyz[m][n][1],1.0/3.0);
      }else{
	xyz[m][n][1]=(pow(29.0/3.0,3)*xyz[m][n][1]+16)/116;
      }
      if(xyz[m][n][2]>0.008856){
	xyz[m][n][2]=pow(xyz[m][n][2],1.0/3.0);
      }else{
	xyz[m][n][2]=(pow(29.0/3.0,3)*xyz[m][n][2]+16)/116;
      }
      //printf("%lf %lf %lf\n",xyz[m][n][0],xyz[m][n][1],xyz[m][n][2]);
      lab[m][n][0]=116*xyz[m][n][1]-16;
      lab[m][n][1]=500*(xyz[m][n][0]-xyz[m][n][1]);
      lab[m][n][2]=200*(xyz[m][n][1]-xyz[m][n][2]);
      //printf("%lf %lf %lf\n",lab[m][n][0],lab[m][n][1],lab[m][n][2]);
    }
  }
}

void convert_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3], double lab3[][WIDTH][3])
{
  int m,n;
  double sum_l=0,sum2_l=0;
  double sum_a=0,sum2_a=0;
  double sum_b=0,sum2_b=0;
  double ave_l,ave2_l;
  double ave_a,ave2_a;
  double ave_b,ave2_b;
  double hensa_l,hensa2_l;
  double hensa_a,hensa2_a;
  double hensa_b,hensa2_b;

  for(m=HEIGHT_START;m<HEIGHT_END;m++){
    for(n=WIDTH_START;n<WIDTH_END;n++){
       sum_l+=lab[m][n][0];
       sum_a+=lab[m][n][1];
       sum_b+=lab[m][n][2];
      }
    }

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
     for(n=WIDTH_START2;n<WIDTH_END2;n++){
       sum2_l+=lab2[m][n][0];
     sum2_a+=lab2[m][n][1];
         sum2_b+=lab2[m][n][2];
   }
    }

    ave_l=sum_l/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
    ave2_l=sum2_l/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));
   ave_a=sum_a/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
   ave2_a=sum2_a/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));
   ave_b=sum_b/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
   ave2_b=sum2_b/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));

  printf("%lf %lf %lf %lf %lf %lf\n",ave_l,ave_a,ave_b,ave2_l,ave2_a,ave2_b);

  sum_l=0;
  sum2_l=0;
  sum_a=0;
  sum2_a=0;
   sum_b=0;
   sum2_b=0;

   for(m=HEIGHT_START;m<HEIGHT_END;m++){
    for(n=WIDTH_START;n<WIDTH_END;n++){
      sum_l+=pow(lab[m][n][0]-ave_l,2);
      sum_a+=pow(lab[m][n][1]-ave_a,2);
     sum_b+=pow(lab[m][n][2]-ave_b,2);
     }
   }

  for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
     for(n=WIDTH_START2;n<WIDTH_END2;n++){
      sum2_l+=pow(lab2[m][n][0]-ave2_l,2);
      sum2_a+=pow(lab2[m][n][1]-ave2_a,2);
     sum2_b+=pow(lab2[m][n][2]-ave2_b,2);
   }
   }

  hensa_l=sqrt(sum_l/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START)));
  hensa2_l=sqrt(sum2_l/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2)));
  hensa_a=sqrt(sum_a/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START)));
  hensa2_a=sqrt(sum2_a/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2)));
  hensa_b=sqrt(sum_b/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START)));
  hensa2_b=sqrt(sum2_b/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2)));

  printf("%lf %lf %lf %lf %lf %lf\n",hensa_l,hensa2_l,hensa_a,hensa2_a,hensa_b,hensa2_b);

  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      //printf("%lf ",lab[m][n][1]);
      lab[m][n][0]-=ave_l;
      //printf("%lf ",ave_a);
      lab[m][n][1]-=ave_a;
      lab[m][n][2]-=ave_b;
      lab[m][n][0]=(hensa2_l/hensa_l)*lab[m][n][0];
      //printf("%lf %lf ",hensa_a,hensa2_a);
      lab[m][n][1]=(hensa2_a/hensa_a)*lab[m][n][1];
      if(lab[m][n][1]<-128){
	lab[m][n][1]=-128;
      }
      if(lab[m][n][1]>128){
	lab[m][n][1]=128;
      }
      //printf("%lf ",lab[m][n][1]);
      lab[m][n][2]=(hensa2_b/hensa_b)*lab[m][n][2];
      if(lab[m][n][2]<-128){
	lab[m][n][2]=-128;
      }
      if(lab[m][n][2]>128){
	lab[m][n][2]=128;
      }
      lab[m][n][0]+=ave2_l;
      if(lab[m][n][0]<0){
	lab[m][n][0]=0;
      }
      if(lab[m][n][0]>100){
	lab[m][n][0]=100;
      }
      lab[m][n][1]+=ave2_a;
      if(lab[m][n][1]<-128){
	lab[m][n][1]=-128;
      }
      //printf("%lf %lf\n",ave2_a,lab[m][n][1]);
      lab[m][n][2]+=ave2_b;
      if(lab[m][n][2]<-128){
	lab[m][n][2]=-128;
      }
      //printf("%lf %lf\n",ave2_a,lab[m][n][1]);
      lab3[m][n][0]=lab[m][n][0];
      lab3[m][n][1]=lab[m][n][1];
      lab3[m][n][2]=lab[m][n][2];
      //printf("%lf %lf %lf\n",lab3[m][n][0],lab3[m][n][1],lab3[m][n][2]);
    }
  }
}

void lab3rgb(double lab[][WIDTH][3], unsigned char rgb[][WIDTH][3])
{
  int m,n;
  double xyz[HEIGHT][WIDTH][3];
  double rgb2[HEIGHT][WIDTH][3];

  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      xyz[m][n][1]=(lab[m][n][0]+16)/116;
      xyz[m][n][0]=(lab[m][n][1]/500)+xyz[m][n][1];
      xyz[m][n][2]=xyz[m][n][1]-(lab[m][n][2]/200);

      if(xyz[m][n][0]>0.206896){
	xyz[m][n][0]=pow(xyz[m][n][0],3.0);
      }else{
	xyz[m][n][0]=pow(3.0/29.0,3)*(116*xyz[m][n][0]-16);
      }
      if(xyz[m][n][1]>0.206896){
	xyz[m][n][1]=pow(xyz[m][n][1],3.0);
      }else{
	xyz[m][n][1]=pow(3.0/29.0,3)*(116*xyz[m][n][1]-16);
      }
      if(xyz[m][n][2]>0.206896){
	xyz[m][n][2]=pow(xyz[m][n][2],3.0);
      }else{
          xyz[m][n][2]=pow(3.0/29.0,3)*(116*xyz[m][n][2]-16);
      }

      rgb2[m][n][0]=3.134187*xyz[m][n][0]-1.617209*xyz[m][n][1]+-0.490694*xyz[m][n][2];
      rgb2[m][n][1]=-0.978749*xyz[m][n][0]+1.916130*xyz[m][n][1]+0.033433*xyz[m][n][2];
      rgb2[m][n][2]=0.071964*xyz[m][n][0]-0.228994*xyz[m][n][1]+1.405754*xyz[m][n][2]; 

      if(rgb2[m][n][0]<=0.003130805){
	rgb2[m][n][0]*=12.92;
      }else{
	rgb2[m][n][0]=pow(1.055*rgb2[m][n][0],1.0/2.4)-0.055;
      }
      if(rgb2[m][n][1]<=0.003130805){
	rgb2[m][n][1]*=12.92;
      }else{
	rgb2[m][n][1]=pow(1.055*rgb2[m][n][1],1.0/2.4)-0.055;
      }
      if(rgb2[m][n][2]<=0.003130805){
	rgb2[m][n][2]*=12.92;
      }else{
	rgb2[m][n][2]=pow(1.055*rgb2[m][n][2],1.0/2.4)-0.055;
      } 

     if(rgb2[m][n][0]<0){
         rgb[m][n][0]=0;
     }else{
         rgb2[m][n][0]=rgb2[m][n][0]*255;
         rgb[m][n][0]=rgb2[m][n][0];
         if(rgb2[m][n][0]>255){
           rgb[m][n][0]=255;
         }
     }
     if(rgb2[m][n][1]<0){
         rgb[m][n][1]=0;
     }else{
         rgb2[m][n][1]=rgb2[m][n][1]*255;
         rgb[m][n][1]=rgb2[m][n][1];
         if(rgb2[m][n][1]>255){
             rgb[m][n][1]=255;
         }
     }
     if(rgb2[m][n][2]<0){
         rgb[m][n][2]=0;
     }else{
         rgb2[m][n][2]=rgb2[m][n][2]*255;
         rgb[m][n][2]=rgb2[m][n][2];
         if(rgb2[m][n][2]>255){
             rgb[m][n][2]=255;
         }
     }
    }
  }
}

void gradation_plate_lab(double lab[][WIDTH][3], double *a, double *b, double *hensa, double *hensa2, int x, int y, int z, int w)
{
    int m,n,i=0,j,count=0,h_count=0,count_a=0,count_b=0,number[(z-x)*(w-y)],number_count=0,tmp;
    double sum_a=0,sum_b=0,sum_c=0;
    double ave_a,ave_b;
    double bunbo=0,bunsi=0;

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
	  sum_a+=pow(lab[m][n][1],2);
	  sum_b+=pow(lab[m][n][2],2);
	  sum_c+=lab[m][n][1]*lab[m][n][2];
            count+=1;
            //printf("%lf %lf %lf\n",H,S,V);
        }
    }

    ave_a=sum_a/count;
    ave_b=sum_b/count;

    //for(m=x;m<z;m++){
    // for(n=y;n<w;n++){
    //     bunsi+=fabs(lab[m][n][1]-ave_a)*fabs(lab[m][n][2]-ave_b);
    //      bunbo+=pow(fabs(lab[m][n][1]-ave_a),2);
    //  }
    // }

    *a=0;
    bunsi=-(sum_a-sum_b)+sqrt(pow(sum_a-sum_b,2)+4*pow(sum_c,2));
    bunbo=2*sum_c;  
  
    *a=bunsi/bunbo;
    *b=ave_b-*a*ave_a;
    //printf("%lf %lf\n",bunsi,bunbo);
    printf("y=%lfx\n",*a);

    
    sum_a=0;
    sum_b=0;
    bunbo=0;
    bunsi=0;
    count=0;
    
    if(*a>30){
      if(*hensa2==0){
	for(m=x;m<z;m++){
	  for(n=y;n<w;n++){
	    sum_a+=lab[m][n][2];
	    sum_b+=lab[m][n][1];
	    count+=1;
	  }
	}

	ave_a=sum_a/count;
	ave_b=sum_b/count;


	bunsi=-(sum_a-sum_b)+sqrt(pow(sum_a-sum_b,2)+4*pow(sum_c,2));
	bunbo=2*sum_c;

	*a=bunsi/bunbo;
	printf("y'=%lfx\n",*a);
	*hensa=1;
      }
    }
}

void gradation_plate_lab2(double lab[][WIDTH][3], double *a, double *b, double *hensa, double *hensa2, int x, int y, int z, int w)
{
    int m,n,i=0,j,count=0,h_count=0,count_a=0,count_b=0,number[(z-x)*(w-y)],number_count=0,tmp;
    double sum_a=0,sum_b=0,sum_c=0;
    double ave_a,ave_b;
    double bunbo=0,bunsi=0;

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
	  sum_a+=pow(lab[m][n][2],2);
	  sum_b+=pow(lab[m][n][1],2);
	  sum_c+=lab[m][n][1]*lab[m][n][2];
            count+=1;
            //printf("%lf %lf %lf\n",H,S,V);
        }
    }

    ave_a=sum_a/count;
    ave_b=sum_b/count;

    //for(m=x;m<z;m++){
    // for(n=y;n<w;n++){
    //     bunsi+=fabs(lab[m][n][1]-ave_a)*fabs(lab[m][n][2]-ave_b);
    //      bunbo+=pow(fabs(lab[m][n][1]-ave_a),2);
    //  }
    // }

    *a=0;
    bunsi=-(sum_a-sum_b)+sqrt(pow(sum_a-sum_b,2)+4*pow(sum_c,2));
    bunbo=2*sum_c;  
  
    *a=bunsi/bunbo;
    *b=ave_b-*a*ave_a;
    //printf("%lf %lf\n",bunsi,bunbo);
    printf("y'=%lfx\n",*a);

}

void gradation_plate_move_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3],  double *a1, double *b1, double *a2, double *b2, double *hensa1, double *hensa2, double *hensa3, double *hensa4, double lab3[][WIDTH][3], double lab4[][WIDTH][3])
{
    int m,n,h,com1,com2,count_a=0,count_b=0;
    double h1,h2;
    double a=0,b=0;
    double katamuki;
    double dista;
    double x,y;
    double dista2=0;
    double dista3=0;
    double g1=0,g2=0;

    for(m=HEIGHT_START;m<HEIGHT_END;m++){
        for(n=WIDTH_START;n<WIDTH_END;n++){
            a+=lab[m][n][1];
            //printf("%lf\t",lab[m][n][1]);
        }
    }

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
        for(n=WIDTH_START2;n<WIDTH_END2;n++){
            b+=lab4[m][n][1];
        }
    }

    a=a/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
    b=b/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));
    //printf("%lf %lf\n",a,b);

    for(m=HEIGHT_START;m<HEIGHT_END;m++){
      for(n=WIDTH_START;n<WIDTH_END;n++){
	dista2+=sqrt(pow(lab[m][n][1],2)+pow(lab[m][n][2],2));
	g1+=sqrt(pow(lab[m][n][0],2)+pow(sqrt(pow(lab[m][n][1],2)+pow(lab[m][n][2],2)),2));
      }
    }

    dista2/=(HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START);
    g1/=(HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START);

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
      for(n=WIDTH_START2;n<WIDTH_END2;n++){
	dista3+=sqrt(pow(lab4[m][n][1],2)+pow(lab4[m][n][2],2));
	g2+=sqrt(pow(lab4[m][n][0],2)+pow(sqrt(pow(lab4[m][n][1],2)+pow(lab4[m][n][2],2)),2));
      }
    }

    dista3/=(HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2);
    g2/=(HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2);
    //printf("%lf %lf\n",dista2,dista3);

    for(m=0;m<HEIGHT;m++){
      for(n=0;n<WIDTH;n++){
            katamuki=(lab[m][n][2]/lab[m][n][1])-*a1+*a2;
            //printf("%lf %lf %lf %lf %lf\n",lab[m][n][2],lab[m][n][1],*a1,*a2,katamuki);
            dista=sqrt(pow(lab[m][n][1],2)+pow(lab[m][n][2],2))*(dista3/dista2);

            x=sqrt(pow(dista,2)/(pow(katamuki,2)+1));

            if(katamuki<0 && (a>0 && b<0) && *a1>0){
	      x*=-1;
            }
	    if(katamuki<0 && (a>0 && b>0) && (*a1>0 && *a2>0)){
	      x*=-1;
            }
	    if(katamuki<0 && (a>0 && b<0) && *a1<0){
	      x*=-1;
            }
	    if(katamuki<0 && (a<0 && b>0) && *a1<0){
	      x*=-1;
            }
            y=katamuki*x;
            lab2[m][n][1]=x;
            lab2[m][n][2]=y;
            lab2[m][n][0]=lab3[m][n][0];
        }
    }
}

void gradation_plate_move_lab2(double lab[][WIDTH][3], double lab2[][WIDTH][3],  double *a1, double *b1, double *a2, double *b2, double *hensa1, double *hensa2, double *hensa3, double *hensa4, double lab3[][WIDTH][3], double lab4[][WIDTH][3])
{
    int m,n,h,com1,com2,count_a=0,count_b=0;
    double h1,h2;
    double a=0,b=0;
    double katamuki;
    double dista;
    double x,y;
    double dista2=0;
    double dista3=0;
    double g1=0,g2=0;

    for(m=HEIGHT_START;m<HEIGHT_END;m++){
        for(n=WIDTH_START;n<WIDTH_END;n++){
            a+=lab[m][n][2];
            //printf("%lf\t",lab[m][n][1]);
        }
    }

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
        for(n=WIDTH_START2;n<WIDTH_END2;n++){
            b+=lab4[m][n][2];
        }
    }

    a=a/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
    b=b/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));
    //printf("%lf %lf\n",a,b);

    for(m=HEIGHT_START;m<HEIGHT_END;m++){
      for(n=WIDTH_START;n<WIDTH_END;n++){
	dista2+=sqrt(pow(lab[m][n][2],2)+pow(lab[m][n][1],2));
	g1+=sqrt(pow(lab[m][n][0],2)+pow(sqrt(pow(lab[m][n][1],2)+pow(lab[m][n][2],2)),2));
      }
    }

    dista2/=(HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START);
    g1/=(HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START);

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
      for(n=WIDTH_START2;n<WIDTH_END2;n++){
	dista3+=sqrt(pow(lab4[m][n][1],2)+pow(lab4[m][n][2],2));
	g2+=sqrt(pow(lab4[m][n][0],2)+pow(sqrt(pow(lab4[m][n][1],2)+pow(lab4[m][n][2],2)),2));
      }
    }

    dista3/=(HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2);
    g2/=(HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2);
    //printf("%lf %lf\n",dista2,dista3);

    for(m=0;m<HEIGHT;m++){
      for(n=0;n<WIDTH;n++){
            katamuki=(lab[m][n][1]/lab[m][n][2])-*a1+*a2;
            //printf("%lf %lf %lf %lf %lf\n",lab[m][n][2],lab[m][n][1],*a1,*a2,katamuki);
            dista=sqrt(pow(lab[m][n][1],2)+pow(lab[m][n][2],2))*(dista3/dista2);

            x=sqrt(pow(dista,2)/(pow(katamuki,2)+1));

            if(katamuki<0 && (a>0 && b<0) && *a1>0){
	      x*=-1;
            }
	    if(katamuki<0 && (a>0 && b>0) && (*a1>0 && *a2>0)){
	      x*=1;
            }
	    if(katamuki<0 && (a>0 && b<0) && *a1<0){
	      x*=-1;
            }
	    if(katamuki<0 && (a<0 && b>0) && *a1<0){
	      x*=-1;
            }
            y=katamuki*x;
            lab2[m][n][1]=y;
            lab2[m][n][2]=x;
            lab2[m][n][0]=lab3[m][n][0];
	    //printf("%lf %lf %lf %lf\n",lab[m][n][2],lab[m][n][1],*a2,katamuki-*a2);
      }
    }
}

