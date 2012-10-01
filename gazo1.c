//4I 4番 伊藤慶二朗
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ppm.h>

extern void color_view(pixel**,int,int,int);

pixel** get_data(char *name,int *col,int *row,gray *mv)
{
  pixel **data;
  FILE *file;

  file=fopen(name,"rb");
  data=ppm_readppm(file,col,row,mv);
  fclose(file);
  return data;
}

double** mat_create()
{
  int y;
  double **mat;

  mat=(double**)malloc(sizeof(double*)*3);

  for(y=0;y<3;y++)
    mat[y]=(double*)calloc(3,sizeof(double));
  
  return mat;
}

void mat_dispose(double **mat)
{
  int y;
  
  for(y=0;y<3;y++)
    free(mat[y]);
  
  free(mat);
}

void mat_init(double **mat)
{
  int x,y;

  for(y=0;y<3;y++)
    {
      for(x=0;x<3;x++)
	mat[y][x]=(x==y ? 1.0 : 0.0);
    }
}

double** mat_copy(double **mat)
{
  int i;
  double **newmat;

  newmat=mat_create();

  for(i=0;i<9;i++)
    newmat[i/3][i%3]=mat[i/3][i%3];

  return newmat;
}

void mat_exch_row(double **mat,int r1,int r2)
{
  int x;
  double temp;

  for(x=0;x<3;x++)
    {
      temp=mat[r1][x];
      mat[r1][x]=mat[r2][x];
      mat[r2][x]=temp;
    }
}

int mat_judge(double **mat)
{
  int mati[3][3],i,n=0;

  for(i=0;i<9;i++)
    mati[i/3][i%3]=(int)(mat[i/3][i%3]*1000);

  for(i=0;i<3;i++)
    {
      n+=mati[0][i]*mati[1][(i+1)%3]*mati[2][(i+2)%3];
      n-=mati[0][2-i]*mati[1][(4-i)%3]*mati[2][(3-i)%3];
    }
  return(n==0 ? 0 : 1);
}

void mat_row_sort_for_g(double **i,double **o)
{
  int x,y;
  
  for(x=0;x<3;x++)
    {
      if((int)(i[x][x]*1000)==0)
	{
	  for(y=0;y<3;y++)
	    {
	      if((int)(i[y][x]*1000)!=0 && !(x>y && (int)(i[x][y]*1000)==0))
		{
		  mat_exch_row(i,x,y);
		  mat_exch_row(o,x,y);
		  break;
		}
	    }
	}
    }
}

double** mat_gauss_j(double **mat)
{
  int px,py,x,y;
  double t,**result;

  result=mat_create();
  mat_init(result);

  for(py=0;py<3;py++)
    {
      mat_row_sort_for_g(mat,result);
      for(y=0;y<3;y++)
	{
	  t=mat[y][y];
	  for(px=3-1;px>=0;px--)
	    {
	      result[y][px]/=t;
	      if(px<py) continue;
	      mat[y][px]/=t;
	    }
	}
      for(y=0;y<3;y++)
	{
	  if(py==y) continue;
	  t=mat[y][py];
	  for(x=3-1;x>=0;x--)
	    {
	      result[y][x]-=result[py][x]*t;
	      if(x<py) continue;
	      mat[y][x]-=mat[py][x]*t;
	    }
	}
    }
  return result;
}

double** mat_reverse(double **mat)
{
  double **rmat,**matcp;

  if(!mat_judge(mat))
    {return (double**)NULL;}
  matcp=mat_copy(mat);
  rmat=mat_gauss_j(matcp);
  mat_dispose(matcp);

  return rmat;
}

void mat_mult(double **mat2,double **mat1)
{
  int x,y,p;
  double **mat1cp;

  mat1cp=mat_copy(mat1);

  for(y=0;y<3;y++)
    for(x=0;x<3;x++)
      {
	mat1[y][x]=0.0;
	for(p=0;p<3;p++)
	  mat1[y][x]+=mat2[y][p]*mat1cp[p][x];
      }
  mat_dispose(mat1cp);
}

void mat_mult_vec(double **mat,double u,double v,double *x,double *y)
{
  int i;
  double n;

  *x=0.0;
  *y=0.0;

  for(i=0;i<3;i++)
    {
      n=(i==0 ? u : (i==1 ? v : 1.0));
      *x+=mat[0][i]*n;
      *y+=mat[1][i]*n;
    }
}

void stretching(double **y)
{
  int a,b;
  double **E;

  printf("横の倍率(％)?\n");
  while(1)
    {
      scanf("%d",&a);
      if(a<10 || a>500)
	printf("10~500(％)で入力してください\n");
      else
	break;
    }
  
  printf("縦の倍率(％)?\n->");
  while(1)
    {
      scanf("%d",&b);
      if(b<10 || b>500)
	printf("10~500(％)で入力してください\n");
      else
	break;
    }

  E=mat_create();
  mat_init(E);

  E[0][0]=a/100.0;
  E[1][1]=b/100.0;

  mat_mult(E,y);
  mat_dispose(E);
}

void rotation(double **y)
{
  int c;
  double **E,sin_c,cos_c;

  printf("角度(deg)?\n");
  while(1)
    {
      scanf("%d",&c);
      if(c<-180 || c>360)
	printf("-180~360(deg)で入力してください\n");
      else
	break;
    }
  E=mat_create();
  mat_init(E);

  sin_c=sin(c*M_PI/180.0);
  cos_c=cos(c*M_PI/180.0);

  E[0][0]=cos_c;
  E[0][1]=sin_c;
  E[1][0]=-sin_c;
  E[1][1]=cos_c;

  mat_mult(E,y);
  mat_dispose(E);
}

void shear(double **y)
{
  int x,d;
  double **E;

  printf("1:水平 2:鉛直\n");
  while(1)
    {
      scanf("%d",&x);
      if(x<1 || x>2)
	printf("1か2で入力してください\n");
      else
	break;
    }

  printf("角度(deg)?\n");
  while(1)
    {
      scanf("%d",&d);
      if(d<=-90 || d>=90)
	printf("-90~90(°)未満で入力してください\n");
      else
	break;
    }

  E=mat_create();
  mat_init(E);

  E[x-1][2-x]=tan(d*M_PI/180.0);

  mat_mult(E,y);
  mat_dispose(E);
}

void shift(double **y)
{
  int e,f;
  double **E;

  printf("横の移動ピクセル?(右が正)\n");
  scanf("%d",&e);

  printf("縦の移動ピクセル?(下が正)\n");
  scanf("%d",&f);
  E=mat_create();
  mat_init(E);

  E[0][2]=e;
  E[1][2]=f;

  mat_mult(E,y);
  mat_dispose(E);
}

void select_transforms(double **y)
{
  int x;
 
  while(1)
    {
      printf("変換処理を選択?　0:決定 1:拡大・縮小 2:回転 3:せん断 4:平行移動\n");
      while(1)
	{
	  scanf("%d",&x);
	  if(x<0 || x>4)
	    printf("0~4で入力してください\n");
	  else
	    break;
	}
      if(x!=0)
	{
	  if(x<3)
	    {
	      if(x==1)
		stretching(y);
	      else
		rotation(y);
	    }
	  else
	    {
	      if(x==3)
		shear(y);
	      else
		shift(y);
	    }
	}
      else
	break;
    }
}

int select_intpol()
{
  int x;

  printf("内挿処理を選択? 1:ニアレストネイバー 2:バイリニア\n");
  while(1)
    {
      scanf("%d",&x);
      if(x<1 || x>2)
	printf("1か2で入力してください\n");
    else
      break;
    }
  return x;
}

void nearest_nb(double **mat,double **rmat,pixel **org_data,int org_col,int org_row,pixel **new_data,int new_col,int new_row)
{
  int u,v,x,y;
  double dx,dy,sx,sy;

  sx=(mat[0][0]<0.0 ? (org_col-1)*mat[0][0] : 0.0)+(mat[0][1]<0.0 ? (org_row-1)*mat[0][1] : 0.0);

  sy=(mat[1][0]<0.0 ? (org_col-1)*mat[1][0] : 0.0)+(mat[1][1]<0.0 ? (org_row-1)*mat[1][1] : 0.0);

  for(v=0;v<new_row;v++)
    {
      for(u=0;u<new_col;u++)
	{
	  mat_mult_vec(rmat,(double)u+sx,(double)v+sy,&dx,&dy);
	  x=(int)(dx+0.5);
	  y=(int)(dy+0.5);
	  if(x<0 || x>=org_col || y<0 || y>=org_row)
	    {
	      new_data[v][u].r=0;
	      new_data[v][u].g=0;
	      new_data[v][u].b=0;
	      continue;
	    }
	  new_data[v][u]=org_data[y][x];
	}
    }
}

void bi_linear(double **mat,double **rmat,pixel **org_data,int org_col,int org_row,pixel **new_data,int new_col,int new_row)
{
  int u,v,x,y,i,lx,uy;
  double dx,dy,sx,sy,r,g,b,rate;
  
  sx=(mat[0][0]<0.0 ? (org_col-1)*mat[0][0] : 0.0)+(mat[0][1]<0.0 ? (org_row-1)*mat[0][1] : 0.0);
  sy=(mat[1][0]<0.0 ? (org_col-1)*mat[1][0] : 0.0)+(mat[1][1]<0.0 ? (org_row-1)*mat[1][1] : 0.0);
  for(v=0;v<new_row;v++)
    {
      for(u=0;u<new_col;u++)
	{
	  mat_mult_vec(rmat,(double)u+sx,(double)v+sy,&dx,&dy);
	  lx=(int)dx+(dx<0.0 && dx-(double)(int)dx ? -1 : 0);
	  uy=(int)dy+(dy<0.0 && dy-(double)(int)dy ? -1 : 0);
	  if(lx+1<0 || lx>=org_col || uy+1<0 || uy>=org_row)
	    {
	      new_data[v][u].r=0;
	      new_data[v][u].g=0;
	      new_data[v][u].b=0;
	      continue;
	    }
	  r=0.0;
	  g=0.0;
	  b=0.0;
	  for(i=0;i<4;i++)
	    {
	      x=lx+i%2;
	      y=uy+i/2;
	      rate=(1.0-fabs((double)x-dx))*(1.0-fabs((double)y-dy));
	      x=(i%2==0 ? (x<0 ? x+1 : x) : (x>=org_col ? x-1 : x));
	      y=(i/2==0 ? (y<0 ? y+1 : y) : (y>=org_row ? y-1 : y));

	      r+=org_data[y][x].r*rate;
	      g+=org_data[y][x].g*rate;
	      b+=org_data[y][x].b*rate;
	    }
	  new_data[v][u].r=(pixval)(r+0.5);
	  new_data[v][u].g=(pixval)(g+0.5);
	  new_data[v][u].b=(pixval)(b+0.5);
	}
    }
}

void image_processing(double **mat,int intpol_type,int org_col,int org_row,gray org_mv,pixel **org_data,int *new_col,int *new_row,gray *new_mv,pixel ***new_data)
{
  int i,j;
  double **rmat;

  *new_col=(int)(fabs((org_col-1)*mat[0][0])+fabs((org_row-1)*mat[0][1])+1.5);
  *new_row=(int)(fabs((org_col-1)*mat[1][0])+fabs((org_row-1)*mat[1][1])+1.5);
  *new_mv=org_mv;
  if(*new_col==0 || *new_row==0)
    {
      printf("エラー 画像サイズが0になりました\n");
      mat_dispose(mat);
      ppm_freearray(org_data,org_row);
      exit(1);
    }
  rmat=mat_reverse(mat);
  if(rmat==(double**)NULL)
    {
      printf("エラー 逆行列の導出に失敗しました\n");
      mat_dispose(mat);
      ppm_freearray(org_data,org_row);
      exit(1);
    }
  *new_data=ppm_allocarray(*new_col,*new_row);
  switch(intpol_type)
    {
    case 1:
      nearest_nb(mat,rmat,org_data,org_col,org_row,*new_data,*new_col,*new_row);
      break;
    case 2:
      bi_linear(mat,rmat,org_data,org_col,org_row,*new_data,*new_col,*new_row);
      break;
    }
  mat_dispose(rmat);
}

int select_output(char *name)
{
  int x;

  printf("出力形式? 1:表示 2:保存 3:保存して表示\n");
  while(1)
    {
      scanf("%d",&x);
      if(x>=1 && x<=3)
	break;
      else
	printf("1~3で入力してください\n");
    }
  if(x!=1)
    {
      printf("保存名?\n");
      scanf("%s",name);
    }
  return x;
}

void output(char *name,int type,pixel **data,int col,int row,gray mv)
{
  FILE *fp;

  if(type==1)
    color_view(data,col,row,(int)mv);
  else
    {
      fp=fopen(name,"wb");
      if(type==2){}
      else
	color_view(data,col,row,(int)mv);
      ppm_writeppm(fp,data,col,row,mv,0);
      printf("保存が完了しました\n");
      fclose(fp);
    }
}

int main(int argc,char *argv[])
{
  int org_col,org_row; gray org_mv; pixel **org_data;
  int new_col,new_row; gray new_mv; pixel **new_data;
  int intpol_type,output_type;
  char o_name[24];
  double **mat;

  if(argc!=2)
    {
      printf("1つの画像名を入力してください\n");
      exit(1);
    }
  org_data=get_data(argv[1],&org_col,&org_row,&org_mv);
  mat=mat_create();
  mat_init(mat);
  select_transforms(mat);
  intpol_type=select_intpol();
  output_type=select_output(o_name);
  image_processing(mat,intpol_type,org_col,org_row,org_mv,org_data,&new_col,&new_row,&new_mv,&new_data);
  mat_dispose(mat);
  ppm_freearray(org_data,org_row);
  output(o_name,output_type,new_data,new_col,new_row,new_mv);
  ppm_freearray(new_data,new_row);

  return 0;
}
