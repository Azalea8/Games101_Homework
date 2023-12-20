//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE, MIRROR, MICROFACET,  Refract, TEST };

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi;} else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        // std::cout << dotProduct(I, N) <<std::endl;
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;

        /*
         B, C其实只需要满足
            B * C = 0
            B * N = 0
            C * N = 0
        */

        // 当 N的 x分量的绝对值大于 N的 y分量的绝对值时，意味着法线 N接近于垂直于 y轴的方向。
        // N的 y分量可能非常小，可能会导致除以接近于 0的数，从而引发数值上的异常。
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            // 这里就用 0来把 N.y消去了
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        // 这里同理
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        // 满足正交矩阵
        B = crossProduct(C, N);

        // 线性代数里面的矩阵相乘
        // 只是框架没有实现矩阵运算，将矩阵运算展开为 标量 乘以 基向量
        return a.x * B + a.y * C + a.z * N;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;

    float roughness;
    float metalness;

    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

    inline float GeometrySchlickGGX(float NdotV, float k);
    inline float DistributionGGX(const Vector3f &N, const Vector3f &H, float a);
    inline float GeometrySmith(float NdotV, float NdotL, float k);
};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1); // [0, 1]
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            // 得到半径为 1一个球的上表面坐标
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);

            // 转换到场景中的坐标系
            return toWorld(localRay, N);
            
            break;
        }
        case MIRROR:
        {
            // 时刻注意方向
            Vector3f localRay = reflect(-wo, N);
            return localRay;
            break;
        }
        case Refract:
        {
            // 时刻注意方向
            Vector3f localRay = refract(-wo, N, ior);
            // std::cout << dotProduct(localRay, N) << std::endl;
            return localRay;
            break;
        }
        case TEST: 
        {
            float r0 = get_random_float();
            float r1 = get_random_float();
            float a2 = roughness * roughness;
            float phi = 2 * M_PI * r1;
            float theta = std::acos(sqrt((1 - r0) / (r0 * (a2 - 1) + 1)));
            // 单位向量半径就直接 1 了，转换为直角坐标系只需要用到 r*sinθ，所以这里直接乘上去了
            float r = std::sin(theta);

            return refract(-wo, toWorld(Vector3f(r * std::cos(phi), r * std::sin(phi), std::cos(theta)), N), ior);
        }
        case MICROFACET:
        {
            
            // 随机一个 ε 和 φ
            float r0 = get_random_float();
            float r1 = get_random_float();
            float a2 = roughness * roughness;
            float phi = 2 * M_PI * r1;
            float theta = std::acos(sqrt((1 - r0) / (r0 * (a2 - 1) + 1)));
            // 单位向量半径就直接 1 了，转换为直角坐标系只需要用到 r*sinθ，所以这里直接乘上去了
            float r = std::sin(theta);
            return reflect(-wo, toWorld(Vector3f(r * std::cos(phi), r * std::sin(phi), std::cos(theta)), N));
        }
    }
}

float Material::pdf(const Vector3f &wo, const Vector3f &wi, const Vector3f &N){

    // 真实的概率密度根本不是这样均匀

    // 有没有发现，我们着色的时候之前要考虑漫反射，镜面反射，环境光
    // 而这里直接让光线在半球面上一顿乱射，也就是我们之前考虑的光线是如何作用于着色点可以看作一个概率问题
    // 镜面反射 (要求 wi wo 关于 N 对称)，这种情况的概率本身就非常低，我们之前还加了个非常大的指数来确保部分高光
    float cosalpha_i_N = dotProduct(wi, N);

    switch(m_type){
        case DIFFUSE: {
            // 也就是为了得到更加真实的渲染，概率密度本身就不该均匀
            // 这个函数的参数附带了 wi，说明它一定是有相关作用的，只是这里的漫反射没有用上

            if (cosalpha_i_N > EPSILON)
                return 0.5 / M_PI;
            else
                return 0.0f;
            break;
        }
        case MIRROR: {
            if (cosalpha_i_N > EPSILON)
                return 1.0f;
            else
                return 0.0f;
            break;
        }
        case Refract: {
            return 1.0f;
                
            break;
        }
        case TEST: 
        {
            if(cosalpha_i_N > EPSILON) {
                // 光线是从物体内部射进空气
                Vector3f h = (ior * wo + wi).normalized(); // h 与 N方向相反
                // std:: cout << dotProduct(h, N) <<std::endl;
                return DistributionGGX(N, -h, roughness) * dotProduct(N, -h) / (4.f * dotProduct(wo, h));
            }
            else if(cosalpha_i_N < -EPSILON) {
                Vector3f h = (ior * wi + wo).normalized();
                return DistributionGGX(N, -h, roughness) * dotProduct(N, -h) / (4.f * dotProduct(wo, -h));
            }
            else {
                return 0.0f;
            }
        }
        case MICROFACET:
        {
            if (cosalpha_i_N > EPSILON) {
                Vector3f h = (wo + wi).normalized();
                return DistributionGGX(N, h, roughness) * dotProduct(N, h) / (4.f * dotProduct(wo, h));
            }else{
                return 0.0f;
            }
            break;
        }
    }
}

inline Vector3f fresnelSchlick(float cosTheta, const Vector3f &F0)
{
    return F0 + (Vector3f(1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float Material::GeometrySchlickGGX(float NdotV, float k)
{
    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}

float Material::GeometrySmith(float NdotV, float NdotL, float k)
{
    float ggx1 = GeometrySchlickGGX(NdotV, k);
    float ggx2 = GeometrySchlickGGX(NdotL, k);

    return ggx1 * ggx2;
}

float Material::DistributionGGX(const Vector3f &N, const Vector3f &H, float a)
{
    float a2     = a*a;
    float NdotH  = std::max(dotProduct(N, H), 0.f);
    float NdotH2 = NdotH*NdotH;
    // std::cout << NdotH2 << std::endl;
    float nom    = a2;
    float denom  = (NdotH2 * (a2 - 1.0) + 1.0);
    // std::cout << denom << std::endl;
    denom        = M_PI * denom * denom;
    // denom = std::max(denom, EPSILON);

    return nom / denom;
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case MIRROR:
        {
            float cosalpha = dotProduct(N, wo);
            float kr;
            if (cosalpha > EPSILON) {
                fresnel(-wi, N, ior, kr);
                Vector3f mirror = 1 / cosalpha;
                return kr * mirror;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case Refract:
        {
            float cosalpha = dotProduct(N, wo);
            float kr;

            fresnel(-wi, N, ior, kr);
            // Vector3f mirror = 1 / std::fabs(cosalpha);
            // std::cout <<  kr * mirror<<std::endl;
            return (1 -kr);
            break;
        }
        case TEST:
        {
            float cosTheta_o_N = dotProduct(N, wo);
            float cosTheta_i_N = dotProduct(N, wi);
            if (cosTheta_i_N * cosTheta_o_N > 0.0f) {
                // std::cout << "}}"<<std::endl;
                // 发生反射
                Vector3f V = wi;
                Vector3f L = wo;
                Vector3f H = (V + L).normalized();
                float NdotV = std::max(dotProduct(N, V), EPSILON);
                float NdotL = std::max(dotProduct(N, L), EPSILON);
                // 直接光照情况下的 k 公式
                float k = (roughness + 1.f) * (roughness + 1.f) / 8.f;
                float D = DistributionGGX(N, H, roughness);
                float G = GeometrySmith(NdotV, NdotL, k);

                Vector3f F0(0.04f);
                F0 = lerp(F0, Kd, metalness);
                Vector3f F = fresnelSchlick(dotProduct(H, V), F0);
                // float F;
                // fresnel(-V, N, ior, F);
                Vector3f fs = D * G * F / (4.f * NdotV  * NdotL);

                // 菲涅尔项就是 ks， kd = 1-ks;
                Vector3f fr =  Kd / M_PI;

                // return (Vector3f(1.0f) - F0) * fr + F0 * fs;
                return (Vector3f(1.0f) - F) * (1 - metalness) * fr + fs;
            }
            else if(cosTheta_i_N * cosTheta_o_N < -EPSILON) {
                // std::cout << "}}"<<std::endl;
                // 发生折射
                Vector3f H;
                float ior_now, ior_i, ior_o;
                if(cosTheta_i_N > 0.0f) {
                    // 光线是从物体内部射进空气
                    H = (ior * wo + wi).normalized(); // h 与 N方向相反
                    ior_now = 1.f / ior;
                    ior_i = 1.f / ior;
                    ior_o = 1;
                }else {
                    H = (ior * wi + wo).normalized(); // h 与 N方向相反
                    ior_now = 1.f / ior;
                    ior_i = 1.f / ior;
                    ior_o = 1;
                }
                
                float NdotWo = std::fabs(dotProduct(N, wo));
                float NdotWi = std::fabs(dotProduct(N, wi));
                float HdotWo = std::fabs(dotProduct(H, wo));
                float HdotWi = std::fabs(dotProduct(H, wi));

                float k = (roughness + 1.f) * (roughness + 1.f) / 8.f;
                float D = DistributionGGX(N, -H, roughness);
                float G = GeometrySmith(NdotWi, NdotWo, k);

                Vector3f F0(0.04f);
                F0 = lerp(F0, Kd, metalness);
                Vector3f F = fresnelSchlick(HdotWi, F0);

                // std::cout << F <<std::endl;
                // float F_;
                // fresnel(-wi, -H, ior, F_);

                // Vector3f F(F_);
                // std::cout << F << std::endl;
                float temp = powf(ior_i * HdotWi + ior_o * HdotWo, 2);

                Vector3f fr = (HdotWo * HdotWi * D * G * (Vector3f(1.0f) - F) * ior_now * ior_now) /  NdotWi * NdotWo;
                // std::cout << D <<std::endl;
                
                Vector3f fs =  Kd / M_PI;;
                return (1- metalness) * fr + F * fs;
                // return (Vector3f(1.0f) - F);
            }
            
            return Vector3f(0.f);
            break;
        }
        case MICROFACET:
        {
            // cosθ是入射光和法线的夹角，也就是光源方向和法线方向的夹角
            float cosTheta = dotProduct(N, wo);
            if(cosTheta > EPSILON) {
                Vector3f V = wi;
                Vector3f L = wo;
                Vector3f H = (V + L).normalized();
                float NdotV = std::max(dotProduct(N, V), EPSILON);
                float NdotL = cosTheta;
                // 直接光照情况下的 k 公式
                float k = (roughness + 1.f) * (roughness + 1.f) / 8.f;
                float D = DistributionGGX(N, H, roughness);
                // std::cout << D <<std::endl;
                float G = GeometrySmith(NdotV, NdotL, k);

                Vector3f F0(0.04f);
                F0 = lerp(F0, Kd, metalness);
                Vector3f F = fresnelSchlick(dotProduct(H, V), F0);
                // float F;
                // fresnel(-V, N, ior, F);
                Vector3f fs = D * G * F / (4.f * NdotV  * NdotL);

                // std::cout << fs <<std::endl;
                // 菲涅尔项就是 ks， kd = 1-ks;
                Vector3f fr =  Kd / M_PI;

                // return (Vector3f(1.0f) - F0) * fr + F0 * fs;
                return (Vector3f(1.0f) - F) * (1 - metalness) * fr + fs;
            }
            return Vector3f(0.f);
            break;
        }
    }
}



#endif //RAYTRACING_MATERIAL_H
