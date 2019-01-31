#ifndef PARPE_COMMON_COST_FUNCTION_H
#define PARPE_COMMON_COST_FUNCTION_H

#include <cassert>
#include <vector>
#include <cmath>

namespace parpe {

/**
 * @brief The CostFunction class.
 *
 * Currently not used. Test implementation to see if we can abstract things a
 * bit more.
 */
class CostFunction {
public:
    virtual ~CostFunction() = default;

    virtual void evaluate(std::vector<double> const& label,
                          std::vector<double> const& prediction,
                          double& cost) const {
        evaluate(label, prediction, 0,
                 std::vector<double * >(0),
                 cost, nullptr);
    }

    virtual void evaluate(std::vector<double> const& label,
                          std::vector<double> const& prediction,
                          int numParameters,
                          std::vector<double *> predictionGradient,
                          double &cost,
                          double *gradient) const = 0;
};

class MeanSquaredError : public CostFunction {
public:
    using CostFunction::evaluate;

    void evaluate(std::vector<double> const& label,
                  std::vector<double> const& prediction,
                  int numParameters,
                  std::vector<double  *> predictionGradient,
                  double &cost,
                  double *gradient) const override {

        assert(label.size() == prediction.size());

        cost = 0.0;

        for(int i = 0; (unsigned) i < label.size(); ++i) {
            cost += std::pow(label[i] - prediction[i], 2);
        }

        cost /= label.size();

        if(gradient) {
            for(int p = 0; p < numParameters; ++p) {
                gradient[p] = 0.0;

                for(int i = 0; (unsigned) i < label.size(); ++i) {
                    gradient[p] += -2.0 * (label[i] - prediction[i])
                            * predictionGradient[i][p];
                }

                gradient[p] /= label.size();
            }
        }
    }
};

}
#endif
