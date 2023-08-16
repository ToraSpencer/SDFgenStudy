#include "makelevelset3.h"
 

// �����߾��룻  find distance x0 is from segment x1-x2
float point_segment_distance(const Vec3f &x0, const Vec3f &x1, const Vec3f &x2)
{
   Vec3f arrow(x2-x1);
   double m2=mag2(arrow);

   // find parameter value of closest point on segment
   float s12=(float)(dot(x2-x0, arrow)/m2);
   if(s12<0) 
      s12=0;
   else if(s12>1) 
      s12=1;
 
   // and find the distance
   return dist(x0, s12*x1+(1-s12)*x2);
}


// �����-����Ƭ���룻   find distance x0 is from triangle x1-x2-x3
float point_triangle_distance(const Vec3f &x0, const Vec3f &x1, const Vec3f &x2, const Vec3f &x3)
{
   // first find barycentric coordinates of closest point on infinite plane
   Vec3f x13(x1-x3), x23(x2-x3), x03(x0-x3);
   float m13=mag2(x13), m23=mag2(x23), d=dot(x13,x23);
   float invdet=1.f/max(m13*m23-d*d,1e-30f);
   float a=dot(x13,x03), b=dot(x23,x03);
   // the barycentric coordinates themselves
   float w23=invdet*(m23*a-d*b);
   float w31=invdet*(m13*b-d*a);
   float w12=1-w23-w31;
   if(w23>=0 && w31>=0 && w12>=0)  // if we're inside the triangle
      return dist(x0, w23*x1+w31*x2+w12*x3); 
    else
    { 
       // we have to clamp to one of the edges
      if(w23>0) // this rules out edge 2-3 for us
         return min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x1,x3));
      else if(w31>0) // this rules out edge 1-3
         return min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x2,x3));
      else // w12 must be >0, ruling out edge 1-2
         return min(point_segment_distance(x0,x1,x3), point_segment_distance(x0,x2,x3));
   }
}


void check_neighbour(const std::vector<Vec3ui> &tris, const std::vector<Vec3f> &vers,
                 SDF_GEN::Array3f &DFvalues, SDF_GEN::Array3i &closest_tris,
                 const Vec3f &ver0, int i0, int j0, int k0, int i1, int j1, int k1)
{
   if(closest_tris(i1,j1,k1)>=0)
   {
      unsigned int vaIdx, vbIdx, vcIdx;
      const Vec3ui& tri = tris[closest_tris(i1, j1, k1)];
      assign(tri, vaIdx, vbIdx, vcIdx);
      float d=point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);

      if(d < DFvalues(i0,j0,k0))
      {
#ifdef USE_MULTITHREADS
          std::lock_guard<std::mutex> guard(g_mutex);
#endif
         DFvalues(i0,j0,k0) = d;
         closest_tris(i0,j0,k0)=closest_tris(i1,j1,k1);
      }
   }
}


void sweep(const std::vector<Vec3ui> &tris, const std::vector<Vec3f> &vers,
    SDF_GEN::Array3f &DFvalues, SDF_GEN::Array3i &closest_tris, const Vec3f &origin, float step,
                  int di, int dj, int dk)
{
    const float ox = origin[0];
    const float oy = origin[1];
    const float oz = origin[2];

    int i0, i1;
   if(di>0)
   {
       i0=1; 
       i1=DFvalues.ni;
   }
   else
   { 
       i0=DFvalues.ni-2; 
       i1=-1; 
   }

   int j0, j1;
   if(dj>0)
   {
       j0=1;
        j1=DFvalues.nj; 
   }
   else
   { 
       j0=DFvalues.nj-2; 
       j1=-1; 
   }

   int k0, k1;
   if(dk>0)
   { 
       k0=1;
       k1=DFvalues.nk;
   }
   else
   { 
       k0=DFvalues.nk-2;
       k1=-1;
   }

   int factor = 0;
   int factorCount = std::ceil(float(k1 - k0)/dk);

   PARALLEL_FOR(0, factorCount, [&](int ft) 
       {
           int k = k0 + dk * ft;
           for (int j = j0; j != j1; j += dj)
           {
               for (int i = i0; i != i1; i += di)
               {
                   Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i - di, j, k);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i, j - dj, k);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i, j, k - dk);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i - di, j, k - dk);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i, j - dj, k - dk);
                   check_neighbour(tris, vers, DFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k - dk);
               }
           }
       });
 
}


// ����(0, 0)(x1, y1)(x2, y2)ȷ���������η������������η���z�����򷵻�1�����򷵻�-1�� ���˻������εĻ�����0;
int orientation(double x1, double y1, double x2, double y2, double &twice_signed_area)
{
   //  ����(0, 0)(x1, y1)(x2, y2)ȷ����������������������������η���z���������������Ϊ������֮Ϊ����
   twice_signed_area=y1*x2-x1*y2;
   if(twice_signed_area>0) 
       return 1;
   else if(twice_signed_area<0) 
       return -1;
   else if(y2>y1) 
       return 1;
   else if(y2<y1) 
       return -1;
   else if(x1>x2) 
       return 1;
   else if(x1<x2) 
       return -1;
   else 
       return 0; // only true when x1==x2 and y1==y2
}


// 2Dƽ���ϲ��Ե�(x0, y0)�Ƿ��� (x1,y1)(x2,y2)(x3,y3)����Χ�ɵ��������ڣ�����ǣ��򷵻�true������abc����������д��õ����������е��������ꣻ
bool point_in_triangle_2d(double x0, double y0, 
                                 double x1, double y1, double x2, double y2, double x3, double y3,
                                 double& a, double& b, double& c)
{
   x1-=x0; x2-=x0; x3-=x0;
   y1-=y0; y2-=y0; y3-=y0;
   int signa=orientation(x2, y2, x3, y3, a);
   if(signa==0) 
       return false;

   int signb=orientation(x3, y3, x1, y1, b);
   if(signb!=signa) 
       return false;
   
   int signc=orientation(x1, y1, x2, y2, c);
   
   if(signc!=signa)
       return false;
   double sum=a+b+c;
   assert(sum!=0);           // if the SOS signs match and are nonkero, there's no way all of a, b, and c are zero.
   a/=sum;
   b/=sum;
   c/=sum;
   return true;
}


void make_level_set3(const std::vector<Vec3ui> &tris, const std::vector<Vec3f> &vers,
                     const Vec3f &startPos, float step, int ni, int nj, int nk,
                    SDF_GEN::Array3f &DFvalues, const int exact_band)
{
    /*
    void make_level_set3(
            const std::vector<Vec3ui> &tris,            ��������Ƭ
            const std::vector<Vec3f> &vers,             ����㼯
            const Vec3f &origin,                                ����դ��ԭ��
            float step,                                                 ����դ�񲽳�
            int ni, int nj, int nk,                                   ����դ������ά�ȵĲ�����
            SDF_GEN::Array3f &DFvalues,                                  ����ľ��볡��ά����
            const int exact_band                                Ĭ��ֵΪ1��������Ƭ��Χ��ʱ������դ������
            )
    
    
    */
    //  ������꣺
    const float ox = startPos[0];
    const float oy = startPos[1];
    const float oz = startPos[2];

   DFvalues.resize(ni, nj, nk);
   DFvalues.assign((ni+nj+nk)*step);                                    // upper bound on distance
   SDF_GEN::Array3i closest_tris(ni, nj, nk, -1);                     // ÿ������դ������������Ƭ����������ʼ��Ϊ-1��
   SDF_GEN::Array3i intersection_count(ni, nj, nk, 0);       // �жϾ������ʱ���ã�intersection_count(i,j,k) is # of tris intersections in (i-1,i]vers{j}vers{k}

   // we begin by initializing distances near the mesh, and figuring out intersection counts

   // 1. ������Ƭ�ı���
    tiktok& tt = tiktok::getInstance();
    tt.start();
 
   PARALLEL_FOR(0, tris.size(), [&](unsigned triIdx)
       {
           unsigned int vaIdx, vbIdx, vcIdx;                  // ��ǰ����Ƭ����������������
           assign(tris[triIdx], vaIdx, vbIdx, vcIdx);

           // 1.1 �ö��㵽����դ�����Ĳ�������ʾ�����λ�ã�
           double vax_sc = ((double)vers[vaIdx][0] - ox) / step, vay_sc = ((double)vers[vaIdx][1] - oy) / step, vaz_sc = ((double)vers[vaIdx][2] - oz) / step;
           double vbx_sc = ((double)vers[vbIdx][0] - ox) / step, vby_sc = ((double)vers[vbIdx][1] - oy) / step, vbz_sc = ((double)vers[vbIdx][2] - oz) / step;
           double vcx_sc = ((double)vers[vcIdx][0] - ox) / step, vcy_sc = ((double)vers[vcIdx][1] - oy) / step, vcz_sc = ((double)vers[vcIdx][2] - oz) / step;

           // ������Ƭ�İ�Χ�У�������դ���е�������ʾ����(i0, j0, k0), (i1, j1, k1)ȷ���ĳ�����ռ䣻
           int i0 = clamp(int(min(vax_sc, vbx_sc, vcx_sc)) - exact_band, 0, ni - 1), i1 = clamp(int(max(vax_sc, vbx_sc, vcx_sc)) + exact_band + 1, 0, ni - 1);
           int j0 = clamp(int(min(vay_sc, vby_sc, vcy_sc)) - exact_band, 0, nj - 1), j1 = clamp(int(max(vay_sc, vby_sc, vcy_sc)) + exact_band + 1, 0, nj - 1);
           int k0 = clamp(int(min(vaz_sc, vbz_sc, vcz_sc)) - exact_band, 0, nk - 1), k1 = clamp(int(max(vaz_sc, vbz_sc, vcz_sc)) + exact_band + 1, 0, nk - 1);

           // ��������Ƭ��Χ���ڵ���������դ��㣬����ǰ����Ƭ�ľ��룬�ҳ���С���룻
           for (int k = k0; k <= k1; ++k) 
               for (int j = j0; j <= j1; ++j) 
                   for (int i = i0; i <= i1; ++i)
           {
               Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);           // ��ǰդ������㣻
               float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);
               if (d < DFvalues(i, j, k))
               {
                   std::lock_guard<std::mutex> guard(g_mutex);
                   DFvalues(i, j, k) = d;                // ����ǰ����ľ����֮ǰ�����С�����¾�������
                   closest_tris(i, j, k) = triIdx;        // 
               }
           }

           // and do intersection counts 
           j0 = clamp((int)std::ceil(min(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);                // ��ǰ����Ƭy��Сֵ��Ӧ��yդ���±ꣻ
           j1 = clamp((int)std::floor(max(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);
           k0 = clamp((int)std::ceil(min(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);            // ��ǰ����Ƭz��Сֵ��Ӧ��zդ���±ꣻ
           k1 = clamp((int)std::floor(max(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);

           for (int k = k0; k <= k1; ++k) 
               for (int j = j0; j <= j1; ++j)
           {
               double x_bc, y_bc, z_bc;       // 
               if (point_in_triangle_2d(j, k, vay_sc, vaz_sc, vby_sc, vbz_sc, vcy_sc, vcz_sc, x_bc, y_bc, z_bc))
               {
                   std::lock_guard<std::mutex> guard(g_mutex);
                   double fi = x_bc * vax_sc + y_bc * vbx_sc + z_bc * vcx_sc;      // intersection i coordinate
                   int i_interval = int(std::ceil(fi));               // intersection is in (i_interval-1,i_interval]

                   if (i_interval < 0)
                       ++intersection_count(0, j, k);          // we enlarge the first interval to include everything to the -vers direction
                   else if (i_interval < ni)
                       ++intersection_count(i_interval, j, k);  // we ignore intersections that are beyond the +vers side of the grid
               }
           }
       });

   tt.endCout("elapsed time of step1: ");

   // 2. and now we fill in the rest of the distances with fast sweeping
   tt.start();
   for(unsigned int pass=0; pass < 2; ++pass)
   {
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, +1, +1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, -1, -1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, +1, -1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, -1, +1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, -1, +1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, +1, -1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, -1, -1);
      sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, +1, +1);
   }
   tt.endCout("elapsed time of step2: ");

   // 3. �����жϣ�then figure out signs (inside/outside) from intersection counts
   tt.start();
   PARALLEL_FOR(0, nk, [&](int k)
       {
           for (int j = 0; j < nj; ++j)
           {
               int total_count = 0;
               // ƽ��x�᷽���ϱ���һ��ֱ���ϵ�դ��ͨ����������ȷ�����볡ֵ�ķ��ţ�
               for (int i = 0; i < ni; ++i)
               {
                   total_count += intersection_count(i, j, k);
                   if (total_count % 2 == 1)                            // if parity of intersections so far is odd,
                   {
                       std::lock_guard<std::mutex> guard(g_mutex);
                       DFvalues(i, j, k) = -DFvalues(i, j, k);            // we are inside the mesh
                   }
               }
           }
       }); 
   tt.endCout("elapsed time of step3: ");
}

