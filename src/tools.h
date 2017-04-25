#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   *
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

private:
    template <class s>
        s cal1 (const s &s1, const s &s2);

    template <class s>
        s cal2(const s &p1,const s &p2, const s &v1 , const s &v2);
};

template <class s>
s Tools::cal1(const s &s1,const s &s2)
{
        return pow(s1,2) + pow(s2,2);
}

template <class s>
s Tools::cal2(const s &p1,const s &p2, const s &v1 , const s &v2)
{
        return p1*(v1*p1 - v2*p2);
}

#endif /* TOOLS_H_ */
