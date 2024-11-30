#pragma once

#include <vector>
#include <cstdio>
#include <cmath>
#include "tracer.h"

enum class OptimizerType {
    Adam,
    BasicGD
};

template<typename ParamType, OptimizerType Type>
class Optimizer {
private:
    std::vector<ParamType> m_gradient; // Gradient
    std::vector<ParamType> m_momentum; // First moment vector
    std::vector<ParamType> m_velocity; // Second moment vector
    ParamType m_alpha; // Learning rate
    ParamType m_beta1; // Decay rate for the first moment estimates
    ParamType m_beta2; // Decay rate for the second moment estimates
    ParamType m_epsilon; // Small value to prevent division by zero
    int step = 0; // Time step
	size_t m_size = 0;

public:
    Optimizer(ParamType alpha, ParamType beta1 = 0.5, ParamType beta2 = 0.999, ParamType epsilon = 1e-8)
        : m_alpha(alpha), m_beta1(beta1), m_beta2(beta2), m_epsilon(epsilon), step(0) {}

	void setSize (size_t num_params_T) {
		m_size = num_params_T;
        m_gradient.clear();
		m_gradient.resize(m_size, ParamType(0));
        if constexpr (Type == OptimizerType::Adam) {
            m_momentum.clear();
            m_momentum.resize(m_size, ParamType(0));
            m_velocity.clear();
            m_velocity.resize(m_size, ParamType(0));
        }
	}

	void set_grad(ParamType val){
		std::fill(m_gradient.begin(), m_gradient.end(), val);
	}

	void zero_grad(){
		std::fill(m_gradient.begin(), m_gradient.end(), ParamType(0));
	}

	std::vector<ParamType> gradient() {
		return m_gradient;
	}

	ParamType* ptr_grad() {
		return &m_gradient[0];
	}

    // Method to update the optimizer with a new gradient
    void update(void* params_ptr) {
        
		ParamType *params= static_cast<ParamType*>(params_ptr);

        // Increment time step
        step++;

        // Update moment vectors and parameters
        for (size_t i = 0; i < m_size; ++i) {
            if constexpr (Type == OptimizerType::Adam) {
                m_momentum[i] = m_beta1 * m_momentum[i] + (1 - m_beta1) * m_gradient[i];
                m_velocity[i] = m_beta2 * m_velocity[i] + (1 - m_beta2) * m_gradient[i] * m_gradient[i];

                ParamType alpha_t = m_alpha * std::sqrt(1 - std::pow(m_beta2, step)) / (1 - std::pow(m_beta1, step));
                params[i] -= alpha_t * m_momentum[i] / (std::sqrt(m_velocity[i]) + m_epsilon);
            } else if constexpr (Type == OptimizerType::BasicGD) {
                params[i] -= m_alpha * m_gradient[i];
            }
        }
    }
};

void render(Cell* opt, float3* out, const RayMarcherExample* pImpl, const float4x4& c2w, const int W, const int H);