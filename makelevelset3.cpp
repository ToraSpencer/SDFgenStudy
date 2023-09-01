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


// �����-����Ƭ���룻   find distance p0 is from triangle p1-p2-p3
float point_triangle_distance(const Vec3f &p0, const Vec3f &p1, const Vec3f &p2, const Vec3f &p3)
{
   // first find barycentric coordinates of closest point on infinite plane
   Vec3f arrow31(p1-p3), arrow32(p2-p3), arrow30(p0-p3);
   float a = dot(arrow31, arrow30);
   float b = dot(arrow32, arrow30);
   float m31=mag2(arrow31), m32=mag2(arrow32), d=dot(arrow31,arrow32);
   float invdet=1.f / max(m31 * m32 - d * d,1e-30f);

   // the barycentric coordinates themselves
   float w32=invdet*(m32*a-d*b);
   float w31=invdet*(m31*b-d*a);
   float w12=1-w32-w31;
   if(w32>=0 && w31>=0 && w12>=0)                // if we're inside the triangle
      return dist(p0, w32*p1+w31*p2+w12*p3); 
    else
    { 
       // we have to clamp to one of the edges
      if(w32>0) // this rules out edge 2-3 for us
         return min(point_segment_distance(p0,p1,p2), point_segment_distance(p0,p1,p3));
      else if(w31>0) // this rules out edge 1-3
         return min(point_segment_distance(p0,p1,p2), point_segment_distance(p0,p2,p3));
      else // w12 must be >0, ruling out edge 1-2
         return min(point_segment_distance(p0,p1,p3), point_segment_distance(p0,p2,p3));
   }
}


void check_neighbour(const std::vector<Vec3ui> &tris, const std::vector<Vec3f> &vers,
                 SDF_GEN::Array3f &SDFvalues, SDF_GEN::Array3i &closest_tris,
                 const Vec3f &ver0, const int i0, const int j0, const int k0, const int i1, const int j1, const int k1)
{
   if(closest_tris(i1, j1, k1) >= 0)
   {
      unsigned int vaIdx, vbIdx, vcIdx;
      const Vec3ui& tri = tris[closest_tris(i1, j1, k1)];
      assign(tri, vaIdx, vbIdx, vcIdx);
      float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);

      if(d < SDFvalues(i0, j0, k0))
      {
#ifdef USE_MULTITHREADS
          std::lock_guard<std::mutex> guard(g_mutex);
#endif
         SDFvalues(i0,j0,k0) = d;
         closest_tris(i0,j0,k0) = closest_tris(i1,j1,k1);
      }
   }
}


void sweep(const std::vector<Vec3ui>&tris, const std::vector<Vec3f> &vers,
    SDF_GEN::Array3f& SDFvalues, SDF_GEN::Array3i &closest_tris, const Vec3f &origin, const float step,
                  const int di, const int dj, const int dk)
{
    const float ox = origin[0];
    const float oy = origin[1];
    const float oz = origin[2];

    int i0, i1, j0, j1, k0, k1;
   if(di > 0)
   {
       i0=1; 
       i1=SDFvalues.ni;
   }
   else
   { 
       i0=SDFvalues.ni-2; 
       i1=-1; 
   }
    
   if(dj > 0)
   {
       j0=1;
        j1=SDFvalues.nj; 
   }
   else
   { 
       j0=SDFvalues.nj-2; 
       j1=-1; 
   }
    
   if(dk > 0)
   { 
       k0=1;
       k1=SDFvalues.nk;
   }
   else
   { 
       k0 = SDFvalues.nk-2;
       k1 = -1;
   }

#ifdef USE_MULTITHREADS
   int factor = 0;
   int factorCount = std::ceil(float(k1 - k0) / dk);

   PARALLEL_FOR(0, factorCount, [&](int ft)
       {
           int k = k0 + dk * ft;
           for (int j = j0; j != j1; j += dj)
           {
               for (int i = i0; i != i1; i += di)
               {
                   Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j, k);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j - dj, k);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j, k - dk);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j, k - dk);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j - dj, k - dk);
                   check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k - dk);
               }
           }
       });
#else
   for (int k = k0; k != k1; k += dk) for (int j = j0; j != j1; j += dj) for (int i = i0; i != i1; i += di)
   {
       Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j, k);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j - dj, k);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j, k - dk);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j, k - dk);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i, j - dj, k - dk);
       check_neighbour(tris, vers, SDFvalues, closest_tris, ver0, i, j, k, i - di, j - dj, k - dk);
   }
#endif
}


// ����(0, 0)(x1, y1)(x2, y2)ȷ���������η������������η���z�����򷵻�1�����򷵻�-1�� ���˻������εĻ�����0;
int orientation(double x1, double y1, double x2, double y2, double &twice_signed_area)
{
   //  ����(0, 0)(x1, y1)(x2, y2)ȷ����������������������������η���z���������������Ϊ������֮Ϊ����
   twice_signed_area=y1*x2-x1*y2;       // norm(arrow1.cross(arrow2))
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


// դ������ϵ�£�2Dƽ���ϲ��Ե�(X0, Y0)�Ƿ��� (X1,Y1), (X2,Y2), (X3,Y3)����Χ�ɵ��������ڣ�����ǣ��򷵻�true������abc����������д��õ����������е��������ꣻ
bool point_in_triangle_2d(double X0, double Y0, 
                                 double X1, double Y1, double X2, double Y2, double X3, double Y3,
                                 double& alpha, double& beta, double& gamma)
{
    // 1. �൱����(X0, Y0)Ϊԭ�㽨���µ�ƽ������ϵ��
   X1-=X0; X2-=X0; X3-=X0;
   Y1-=Y0; Y2-=Y0; Y3-=Y0;

   // 2. ͨ���������������������������(X0, Y0)������Χ�ɵ��������е��������ꣻ
   double S23, S31, S12;
   int signa=orientation(X2, Y2, X3, Y3, S23);
   if(signa==0) 
       return false;

   int signb=orientation(X3, Y3, X1, Y1, S31);
   if(signb!=signa) 
       return false;
   
   int signc=orientation(X1, Y1, X2, Y2, S12);
   
   if(signc!=signa)
       return false;
   double sum = S23 + S31 + S12;
   assert(sum!=0);           // if the SOS signs match and are nonkero, there's no waY all of alpha, beta, and gamma are zero.
   alpha = S23 / sum;
   beta = S31 / sum;
   gamma = S12 / sum;
   return true;
}


void make_level_set3(const std::vector<Vec3ui> &tris, const std::vector<Vec3f> &vers,
                     const Vec3f &startPos, float step, int ni, int nj, int nk,
                    SDF_GEN::Array3f &SDFvalues, const int exact_band)
{
    /*
    void make_level_set3(
            const std::vector<Vec3ui> &tris,            ��������Ƭ
            const std::vector<Vec3f> &vers,             ����㼯
            const Vec3f &origin,                                ����դ��ԭ��
            float step,                                                 ����դ�񲽳�
            int ni, int nj, int nk,                                   ����դ������ά�ȵĲ�����
            SDF_GEN::Array3f &SDFvalues,                                  ����ľ��볡��ά����
            const int exact_band                                Ĭ��ֵΪ1��������Ƭ��Χ��ʱ������դ������
            )
    
    
    */
    //  ������꣺
    const float ox = startPos[0];
    const float oy = startPos[1];
    const float oz = startPos[2];

   SDFvalues.resize(ni, nj, nk);
   SDFvalues.assign((ni+nj+nk)*step);                                    // upper bound on distance
   SDF_GEN::Array3i closest_tris(ni, nj, nk, -1);                     // ÿ������դ������������Ƭ����������ʼ��Ϊ-1��
   SDF_GEN::Array3i intersection_count(ni, nj, nk, 0);       // �жϾ������ʱ���ã�intersection_count(i,j,k) is # of tris intersections in (i-1,i]vers{j}vers{k}

   // we begin by initializing distances near the mesh, and figuring out intersection counts

   // 1. ������Ƭ�ı���
#ifdef USE_MULTITHREADS

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
           for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j) for (int i = i0; i <= i1; ++i)
           {
               Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);           // ��ǰդ������㣻
               float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);
               if (d < SDFvalues(i, j, k))
               {
                   std::lock_guard<std::mutex> guard(g_mutex);
                   SDFvalues(i, j, k) = d;                // ����ǰ����ľ����֮ǰ�����С�����¾�������
                   closest_tris(i, j, k) = triIdx;        // 
               }
           }

           // and do intersection counts 
           j0 = clamp((int)std::ceil(min(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);                // ��ǰ����Ƭy��Сֵ��Ӧ��yդ���±ꣻ
           j1 = clamp((int)std::floor(max(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);
           k0 = clamp((int)std::ceil(min(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);            // ��ǰ����Ƭz��Сֵ��Ӧ��zդ���±ꣻ
           k1 = clamp((int)std::floor(max(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);

           for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j)
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


#else

   for (unsigned int triIdx = 0; triIdx < tris.size(); ++triIdx)
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
       for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j) for (int i = i0; i <= i1; ++i)
       {
           Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);           // ��ǰդ������㣻
           float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);
           if (d < SDFvalues(i, j, k))
           {
               SDFvalues(i, j, k) = d;                // ����ǰ����ľ����֮ǰ�����С�����¾�������
               closest_tris(i, j, k) = triIdx;        // 
           }
       }

       // and do intersection counts
       // int i0 = clamp(int(min(vax_sc, vbx_sc, vcx_sc)) - exact_band, 0, ni - 1)
       j0 = clamp((int)std::ceil(min(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);                // ��ǰ����Ƭy��Сֵ��Ӧ��yդ���±ꣻ
       j1 = clamp((int)std::floor(max(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);
       k0 = clamp((int)std::ceil(min(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);            // ��ǰ����Ƭz��Сֵ��Ӧ��zդ���±ꣻ
       k1 = clamp((int)std::floor(max(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);

       for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j)
       {
           double x_bc, y_bc, z_bc;       // 
           if (point_in_triangle_2d(j, k, vay_sc, vaz_sc, vby_sc, vbz_sc, vcy_sc, vcz_sc, x_bc, y_bc, z_bc))
           {
               double fi = x_bc * vax_sc + y_bc * vbx_sc + z_bc * vcx_sc;      // intersection i coordinate
               int i_interval = int(std::ceil(fi));               // intersection is in (i_interval-1,i_interval]

               if (i_interval < 0)
                   ++intersection_count(0, j, k);          // we enlarge the first interval to include everything to the -vers direction
               else if (i_interval < ni)
                   ++intersection_count(i_interval, j, k);  // we ignore intersections that are beyond the +vers side of the grid
           }
       }
   }
#endif

   // 2. and now we fill in the rest of the distances with fast sweeping 
   for(unsigned int pass=0; pass < 2; ++pass)
   {
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, +1, +1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, -1, -1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, +1, -1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, -1, +1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, -1, +1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, +1, -1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, -1, -1);
      sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, +1, +1);
   } 

   // 3. �����жϣ�then figure out signs (inside/outside) from intersection counts 
#ifdef USE_MULTITHREADS
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
                       SDFvalues(i, j, k) = -SDFvalues(i, j, k);            // we are inside the mesh
                   }
               }
           }
       });
#else
   for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j)
   {
       int total_count = 0;
       // ƽ��x�᷽���ϱ���һ��ֱ���ϵ�դ��ͨ����������ȷ�����볡ֵ�ķ��ţ�
       for (int i = 0; i < ni; ++i)
       {
           total_count += intersection_count(i, j, k);
           if (total_count % 2 == 1)                            // if parity of intersections so far is odd,
               SDFvalues(i, j, k) = -SDFvalues(i, j, k);            // we are inside the mesh
       }
   }
#endif
}

