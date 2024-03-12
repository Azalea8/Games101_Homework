#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    // start -> 起始节点的坐标，end -> 最终节点的坐标，num_nodes -> 节点个数
    // node_mass -> 节点的质量，k -> 弹簧系数
    // pinned_nodes -> 弹簧两个端点的索引
    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        //        Comment-in this part when you implement the constructor
        //        for (auto &i : pinned_nodes) {
        //            masses[i]->pinned = true;
        //        }
        for (int i = 0; i < num_nodes; ++i)
        {
            // 转成double再做除法，否则就会截断成int
            Vector2D pos = start + (end - start) * ((double)i / ((double)num_nodes - 1.0));
            masses.push_back(new Mass(pos, node_mass, false)); // 放质点
        }
        for (int i = 0; i < num_nodes - 1; ++i)
            springs.push_back(new Spring(masses[i], masses[i + 1], k)); // 放绳子
        for (auto &i : pinned_nodes)
            masses[i]->pinned = true;
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node

            // 实现胡克定律：遍历所有的弹簧，根据拉伸的长度计算拉力
            auto len = (s->m1->position - s->m2->position).norm();
            s->m1->forces += -s->k * (s->m1->position - s->m2->position) / len * (len - s->rest_length);
            s->m2->forces += -s->k * (s->m2->position - s->m1->position) / len * (len - s->rest_length);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position

                // 显示欧拉法：下一个时刻位置用当前速度计算
                auto a = m->forces / m->mass + gravity;
                m->position += m->velocity * delta_t; // 先更新位置，就是用上一时刻的速度更新
                m->velocity += a * delta_t;

                /* // 隐式欧拉法：下一个时刻位置用下一时刻速度计算
                auto a = m->forces / m->mass + gravity;
                m->velocity += a * delta_t;
                m->position += m->velocity * delta_t; // 后更新位置，就是用这一时刻的速度更新 */
                // TODO (Part 2): Add global damping
            }
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            auto len = (s->m1->position - s->m2->position).norm();
            s->m1->forces += -s->k * (s->m1->position - s->m2->position) / len * (len - s->rest_length);
            s->m2->forces += -s->k * (s->m2->position - s->m1->position) / len * (len - s->rest_length);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // 显示Verlet法：把弹簧的长度保持原长作为约束，移动每个质点的位置
                Vector2D temp_position = m->position;
                auto a = m->forces / m->mass + gravity;
                // TODO (Part 3.1): Set the new position of the rope mass
                m->position = temp_position + (temp_position - m->last_position) + a * delta_t * delta_t;
                m->last_position = temp_position;
            }
            m->forces = Vector2D(0, 0);
        }
    }
}
