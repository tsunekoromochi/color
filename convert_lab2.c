#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEIGHT        129
#define WIDTH         208
#define HEIGHT_START  50
#define WIDTH_START   70
#define HEIGHT_END    90
#define WIDTH_END     190
#define HEIGHT_START2 40
#define WIDTH_START2  5
#define HEIGHT_END2   85
#define WIDTH_END2    55


void ipr_load_ppm(unsigned char image[][WIDTH][3], const char path[]);
void ipr_save_ppm(unsigned char image[][WIDTH][3], const char path[]);
void ipr_rgb3hsv(unsigned char rgb[][WIDTH][3], int hsv[][WIDTH][3]);
void ipr_hsv3rgb(int hsv[][WIDTH][3], unsigned char rgb[][WIDTH][3]);
void gradation_plate(int hsv[][WIDTH][3], double*, double*, double*, double*, int x, int y, int z, int w);
void gradation_plate_move(int hsv[][WIDTH][3], int plate[][WIDTH][3], double*, double*, double*, double*, double*, double*, double*, double*);
void gradation_plate_move2(int hsv[][WIDTH][3], int plate[][WIDTH][3], double*, double*, double*, double*, double*, double*, double*, double*);
void gradation_plate_juusin(int hsv[][WIDTH][3], double*);
void gradation_plate_convert(int plate[][WIDTH][3], double*, double*);
void gradation_plate_range(unsigned char rgb[][WIDTH][3], unsigned char rgb2[][WIDTH][3], int x, int y, int z, int w);
void rgb3lab(unsigned char rgb[][WIDTH][3], double lab[][WIDTH][3]);
void convert_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3], double lab3[][WIDTH][3]);
void lab3rgb(double lab[][WIDTH][3], unsigned char rgb[][WIDTH][3]);
void gradation_plate_lab(double lab[][WIDTH][3], double*, double*, double*, double*, int x, int y, int z, int w);
void gradation_plate_move_lab(double lab[][WIDTH][3], double lab2[][WIDTH][3], double*, double*, double*, double*, double*, double*, double*, double*, double lab3[][WIDTH][3], double lab4[][WIDTH][3]);

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
    double hensa3;
    double hensa4;
    double hensa5;
    double hensa6;
    double lab2[HEIGHT][WIDTH][3];
    unsigned char rgb_lab2[HEIGHT][WIDTH][3];
    FILE *fp3;
    double lab_plate[HEIGHT][WIDTH][3];
    double lab_plate2[HEIGHT][WIDTH][3];

    if (argc != 8) {
      fprintf(stderr, "Usage: %s source destination\n", argv[0]);
        exit(1);
    }

    ipr_load_ppm(src_img, argv[1]);
    ipr_load_ppm(src_img2, argv[2]);

    ipr_rgb3hsv(src_img, dst_img);
    ipr_rgb3hsv(src_img2, dst_img2);

    printf("対象画像\n");
    gradation_plate(dst_img, &a1, &b1, &c1, &hensa1, HEIGHT_START, WIDTH_START, HEIGHT_END, WIDTH_END);

    printf("参照画像\n");
    gradation_plate(dst_img2, &a2, &b2, &c2, &hensa2, HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);
    
    gradation_plate_move(dst_img, plate, &a1, &b1, &c1, &a2, &b2, &c2, &hensa1, &hensa2);
    gradation_plate_juusin(dst_img, &R1);
    gradation_plate_juusin(dst_img2, &R2);
    gradation_plate_convert(plate, &R1, &R2);

    ipr_hsv3rgb(plate, img);

    gradation_plate_range(src_img, range, HEIGHT_START, WIDTH_START, HEIGHT_END, WIDTH_END);
    gradation_plate_range(src_img2, range2,  HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);

    ipr_save_ppm(img, argv[3]);
    ipr_save_ppm(range, argv[4]);
    ipr_save_ppm(range2, argv[5]);

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

    ipr_save_ppm(rgb_lab, argv[6]);

    gradation_plate_lab(lab_plate, &a3, &b3, &hensa3, &hensa4, HEIGHT_START, WIDTH_START, HEIGHT_END, WIDTH_END);
    gradation_plate_lab(lab_plate2, &a4, &b4, &hensa5, &hensa6,  HEIGHT_START2, WIDTH_START2, HEIGHT_END2, WIDTH_END2);
    gradation_plate_move_lab(lab_plate, lab2, &a3, &b3, &a4, &b4, &hensa3, &hensa5,  &hensa4, &hensa6, img_lab3, lab_plate2);

    lab3rgb(lab2, rgb_lab2);

    fp3=fopen("lab3.txt","w");
    for(m=HEIGHT_START;m<HEIGHT_END;m++){
        for(n=WIDTH_START;n<WIDTH_END;n++){
            fprintf(fp,"%lf %lf %lf\n",lab2[m][n][1],lab2[m][n][2],lab2[m][n][0]);
      }
    }
    fclose(fp3);

    for(m=0;m<1;m++){
      for(n=0;n<WIDTH;n++){
          //printf("%d %d %d\n",rgb_lab[m][n][0],rgb_lab[m][n][1],rgb_lab[m][n][2]);
          //printf("%lf %lf %lf\n",img_lab3[m][n][0],img_lab3[m][n][1],img_lab3[m][n][2]);
          //printf("%lf %lf %lf\n",lab2[m][n][0],lab2[m][n][1],lab2[m][n][2]);
      }
    }

    ipr_save_ppm(rgb_lab2, argv[7]);
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

void ipr_rgb3hsv(unsigned char rgb[][WIDTH][3], int hsv[][WIDTH][3])
{
    int m,n;
    double max,min,H,S,V;
    for(m=0;m<HEIGHT;m++){
        for(n=0;n<WIDTH;n++){
            if(rgb[m][n][0]==rgb[m][n][1] && rgb[m][n][0]==rgb[m][n][2]){
                hsv[m][n][0]=H=0;
                hsv[m][n][1]=S=0;
                hsv[m][n][2]=V=max=min=rgb[m][n][0];
            }else if(rgb[m][n][0]<=rgb[m][n][2] && rgb[m][n][1]<=rgb[m][n][2]){
                if(rgb[m][n][0]<=rgb[m][n][1]){
                    min=rgb[m][n][0];
                }else{
                    min=rgb[m][n][1];
                }
                max=V=rgb[m][n][2];
                S=255*(max-min)/max;
                H=60*((rgb[m][n][0]-rgb[m][n][1])/(max-min))+240;
                if(H<0){
                    H+=360;
                }
                hsv[m][n][0]=H;
                hsv[m][n][1]=S;
                hsv[m][n][2]=V;
            }else if(rgb[m][n][0]<=rgb[m][n][1] && rgb[m][n][2]<=rgb[m][n][1]){
                if(rgb[m][n][0]<=rgb[m][n][2]){
                    min=rgb[m][n][0];
                }else{
                    min=rgb[m][n][2];
                }
                max=V=rgb[m][n][1];
                S=255*(max-min)/max;
                H=60*((rgb[m][n][2]-rgb[m][n][0])/(max-min))+120;
                if(H<0){
                    H+=360;
                }
                hsv[m][n][0]=H;
                hsv[m][n][1]=S;
                hsv[m][n][2]=V;
            }else if(rgb[m][n][1]<=rgb[m][n][0] && rgb[m][n][2]<=rgb[m][n][0]){
                if(rgb[m][n][1]<=rgb[m][n][2]){
                    min=rgb[m][n][1];
                }else{
                    min=rgb[m][n][2];
                }
                max=V=rgb[m][n][0];
                S=255*(max-min)/max;
                H=60*((rgb[m][n][1]-rgb[m][n][2])/(max-min));
                if(H<0){
                    H+=360;
                }
                hsv[m][n][0]=H;
                hsv[m][n][1]=S;             
                hsv[m][n][2]=V;
            }
            //printf("%d %d %d\n",hsv[m][n][0],hsv[m][n][1],hsv[m][n][2]);
        }
    }
}

void ipr_hsv3rgb(int hsv[][WIDTH][3], unsigned char rgb[][WIDTH][3])
{
    int m,n;
    double max,min,H,S,V,R,G,B;
    for(m=0;m<HEIGHT;m++){
      for(n=0;n<WIDTH;n++){
            H=hsv[m][n][0];
            S=hsv[m][n][1];
            V=hsv[m][n][2];
            max=V;
            min=max-((S/255)*max);
            if(H>300){
                R=max;
                G=min;
                B=((360-H)/60)*(max-min)+min;
            }else if(H>240){
                R=((H-240)/60)*(max-min)+min;
                G=min;
                B=max;
            }else if(H>180){
                R=min;
                G=((240-H)/60)*(max-min)+min;
                B=max;
            }else if(H>120){
                R=min;
                G=max;
                B=((H-120)/60)*(max-min)+min;
            }else if(H>60){
                R=((120-H)/60)*(max-min)+min;
                G=max;
                B=min;
            }else{
                R=max;
                G=(H/60)*(max-min)+min;
                B=min;
            }
            rgb[m][n][0]=R;
            rgb[m][n][1]=G;             
            rgb[m][n][2]=B;
            //printf("%d %d %d\n",rgb[m][n][0],rgb[m][n][1],rgb[m][n][2]);
        }
    }
}

void gradation_plate(int hsv[][WIDTH][3], double *a, double *b, double *c, double *hensa, int x, int y, int z, int w)
{
    int m,n,i=0,j,count=0,h_count=0,count_a=0,count_b=0,number[(z-x)*(w-y)],number_count=0,tmp;
    double H,S,V;
    double A=0,B=0,C=0,D=0,E=0,F=0,G=0,I=0;

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
            if(hsv[m][n][0]<10){
                count_a++;
            }
            if(hsv[m][n][0]>350){
                count_b++;
            }
        }
    }

    //printf("%d %d",count_a,count_b);

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
            if(count_a>count_b){
                if(hsv[m][n][0]>350){
                    hsv[m][n][0]-=360;
                }
            }else{
                if(hsv[m][n][0]<10){
                    hsv[m][n][0]+=360;
                }
            }

            H=hsv[m][n][0];
            //printf("%lf\n",H);
            S=hsv[m][n][1];
            V=hsv[m][n][2];
            A+=S*S;
            B+=V*V;
            C+=S*V;
            D+=S*H;
            E+=V*H;
            F+=S;
            G+=V;
            I+=H;
            count+=1;
            //printf("%lf %lf %lf\n",H,S,V);
        }
    }
    *a=-((C*(G*I-count*E)+F*(E*G-B*I)+D*(count*B-G*G))/(A*(G*G-count*B)-2*C*F*G+B*F*F+count*C*C));
    *b=(A*(G*I-count*E)-F*(C*I+D*G)+E*F*F+count*C*D)/(A*(G*G-count*B)-2*C*F*G+B*F*F+count*C*C);
    *c=(A*(E*G-B*I)+C*C*I-C*D*G+(B*D-C*E)*F)/(A*(G*G-count*B)-2*C*F*G+B*F*F+count*C*C);

    printf("%lf %lf %lf\n",*a,*b,*c);

    *hensa=0;
    for(m=x;m<z;m++){
      for(n=y;n<w;n++){
	*hensa+=pow((int)(hsv[m][n][0]-(*a*hsv[m][n][1]+*b*hsv[m][n][2]+*c)),2);
	h_count+=1;
      }
    }
    *hensa/=h_count;
    *hensa=sqrt(*hensa);

}

void gradation_plate_move(int hsv[][WIDTH][3], int plate[][WIDTH][3],  double *a1, double *b1, double *c1, double *a2, double *b2, double *c2, double *hensa1, double *hensa2)
{
    int m,n,h,com1,com2,x,count_a=0,count_b=0;
  double h1,h2;
 
  for(m=0;m<HEIGHT;m++){
      for(n=0;n<WIDTH;n++){
          if(hsv[m][n][0]<10){
              count_a++;
          }
          if(hsv[m][n][0]>350){
              count_b++;
          }
      }
  }

  for(m=0;m<HEIGHT;m++){
      for(n=0;n<WIDTH;n++){
          if(count_a>count_b){
              if(hsv[m][n][0]>350){
                  hsv[m][n][0]-=360;
              }
          }else{
              if(hsv[m][n][0]<10){
                  hsv[m][n][0]+=360;
              }
          }
          h1=*a1*hsv[m][n][1]+*b1*hsv[m][n][2]+*c1;
      h2=*a2*hsv[m][n][1]+*b2*hsv[m][n][2]+*c2;
      com1=abs((int)(hsv[m][n][0]-h1));
      if(hsv[m][n][0]<=h1){
	com2=abs((int)(360+hsv[m][n][0]-h1));
      }else{
	com2=abs((int)(hsv[m][n][0]-h1+360));
      }
      if(com1<=com2){
	x=0;
      }else{
	x=1;
      }
      if(x==0){
	h=(*hensa2/(*hensa1))*(int)(hsv[m][n][0]-h1)+h2;
      }else{
	if(hsv[m][n][0]<=h1){
	  h=(*hensa2/(*hensa1))*(int)(360+hsv[m][n][0]-h1)+h2;
      	}
      }
      if(h<0){
	h=h*(-1);
      }
      plate[m][n][0]=h;
      if(plate[m][n][0]>360){
	plate[m][n][0]=plate[m][n][0]%360;
      } 
      plate[m][n][1]=hsv[m][n][1];
      plate[m][n][2]=hsv[m][n][2];
    }
  }
}

void gradation_plate_juusin(int hsv[][WIDTH][3], double *R)
{
  int m,n,R_count=0;

  *R=0;
  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      if(abs(hsv[m][n][1]-hsv[m][n][2])<0.05){
	*R+=sqrt(pow(hsv[m][n][1],2)+pow(hsv[m][n][1],2));
	R_count+=1;
      }
    }
  }

  *R/=R_count;
}

void gradation_plate_convert(int plate[][WIDTH][3], double *R1, double *R2)
{
  int m,n,tmp;

  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      if(plate[m][n][1]<plate[m][n][2]){
	tmp=plate[m][n][1];
	plate[m][n][1]=(*R2/(*R1))*plate[m][n][1];
	plate[m][n][2]=plate[m][n][2]-tmp+plate[m][n][1];
	if(plate[m][n][1]>255){
	  plate[m][n][1]=255;
	}
	if(plate[m][n][2]>255){
	  plate[m][n][2]=255;
	}
      }else{
	tmp=plate[m][n][2];
	plate[m][n][2]=(*R2/(*R1))*plate[m][n][2];
	plate[m][n][1]=plate[m][n][1]-tmp+plate[m][n][2];
	if(plate[m][n][1]>255){
	  plate[m][n][1]=255;
	}
	if(plate[m][n][2]>255){
	  plate[m][n][2]=255;
	}
      }
    }
  }
}

void gradation_plate_move2(int hsv[][WIDTH][3], int plate[][WIDTH][3],  double *a1, double *b1, double *c1, double *a2, double *b2, double *c2, double *hensa1, double *hensa2)
{
  int m,n,h,com1,com2,x;
  double h1,h2;
 
  for(m=0;m<HEIGHT;m++){
    for(n=0;n<WIDTH;n++){
      h1=*a1*hsv[m][n][1]+*b1*hsv[m][n][2]+*c1;
      if(h1>360){
	h1=(int)h1%360;
      }
      if(abs(hsv[m][n][0]-h1)>320){
	if(hsv[m][n][0]<h1){
	  hsv[m][n][0]+=360;
	}else{
	  h1+=360;
	}
      }
      h2=*a2*plate[m][n][1]+*b2*plate[m][n][2]+*c2;
      if(h2>360){
	h2=(int)h2%360;
      }
      if(abs(hsv[m][n][0]-h2)>320){
	if(hsv[m][n][0]<h2){
	  hsv[m][n][0]+=360;
	}else{
	  h2+=360;
	}
      }
      h=hsv[m][n][0]+h2-h1;
      if(h<0){
	h+=360;
      }
      plate[m][n][0]=h;
      if(plate[m][n][0]>360){
	plate[m][n][0]=plate[m][n][0]%360;
      }
      //printf("%d %d %lf %lf\n",plate[m][n][0],hsv[m][n][0],h2,h1);
     //printf("%lf\n",h2);
    }
  }
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

  //printf("%lf %lf %lf %lf %lf %lf\n",ave_l,ave_a,ave_b,ave2_l,ave2_a,ave2_b);

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
    double sum_a=0,sum_b=0;
    double ave_a,ave_b;
    double bunbo=0,bunsi=0;

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
            sum_a+=lab[m][n][1];
            sum_b+=lab[m][n][2];
            count+=1;
            //printf("%lf %lf %lf\n",H,S,V);
        }
    }

    ave_a=sum_a/count;
    ave_b=sum_b/count;

    //for(m=x;m<z;m++){
    //  for(n=y;n<w;n++){
    //      bunsi+=fabs(lab[m][n][1]-ave_a)*fabs(lab[m][n][2]-ave_b);
    //      bunbo+=pow(fabs(lab[m][n][1]-ave_a),2);
    //  }
    //}

    for(m=x;m<z;m++){
        for(n=y;n<w;n++){
            bunsi+=lab[m][n][1]*lab[m][n][2];
            bunbo+=pow(lab[m][n][1],2);
        }
    }
  
    printf("%lf %lf\n",bunsi,bunbo);
    *a=bunsi/bunbo;
    *b=ave_b-*a*ave_a;

    printf("y=%lfx\n",*a);

    *hensa=0;
    for(m=x;m<z;m++){
      for(n=y;n<w;n++){
          *hensa+=pow((int)(lab[m][n][2]-(*a*lab[m][n][1])),2);
          *hensa2+=pow((int)(lab[m][n][1]-((lab[m][n][2]-*b)/(*a))),2);
	h_count+=1;
      }
    }
    *hensa/=h_count;
    *hensa=sqrt(*hensa);
    *hensa2/=h_count;
    *hensa2=sqrt(*hensa2);

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
        }
    }

    for(m=HEIGHT_START2;m<HEIGHT_END2;m++){
      for(n=WIDTH_START2;n<WIDTH_END2;n++){
	b+=lab4[m][n][1];
      }
    }


    a=a/((HEIGHT_END-HEIGHT_START)*(WIDTH_END-WIDTH_START));
    b=b/((HEIGHT_END2-HEIGHT_START2)*(WIDTH_END2-WIDTH_START2));

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
            y=katamuki*x;
            lab2[m][n][1]=x;
            lab2[m][n][2]=y;
            lab2[m][n][0]=lab[m][n][0];
        }
    }
}

