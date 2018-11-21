#include <iostream>
#include "fusion.h"
#include "gflags/gflags.h"
#include "INIReader.h"
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <ctime>
#include <ratio>
#include <chrono>
#include <fstream>

using namespace nsek::fusion;
using namespace nnty;
using namespace std::chrono;

DEFINE_bool(verbose, false, "Display program name before message");
DEFINE_string(config, "/home/chern/", "Message to print");

struct TaskLogger {

  std::ofstream outFile;

  TaskLogger(const std::string &id) {
    outFile.open(id.c_str(), std::ofstream::out);
  }

  //M->setLogHandler([=](const std::string &msg) { std::cout << msg << std::flush; });

  void log(const std::string &msg) {
    outFile << msg;
  }

  ~TaskLogger() {
    outFile.close();
  }
};

static bool IsNonEmptyMessage(const char *flagname, const std::string &value) {
  return value[0] != '\0';
}

void loadWeight(const std::string &filename, size_t length, std::vector<double> &output) {
  std::ifstream y;
  y.open(filename);
  if (y.is_open()) {
    std::string line;
    while (std::getline(y, line)) {
      output.push_back(boost::lexical_cast<double>(line));
    }
    y.close();
  }

  if (length != output.size()) {
    throw std::runtime_error("length of w is not correct! Abort!");
  }

}

void loady(const std::string &filename, size_t length, std::vector<double> &output) {
  std::ifstream y;
  y.open(filename);
  if (y.is_open()) {
    std::string line;
    while (std::getline(y, line)) {
      output.push_back(boost::lexical_cast<double>(line));
    }
    y.close();
  }

  if (length != output.size()) {
    std::ostringstream oss;
    oss << "length of y is not correct! Abort! expecting " << length << ", actually " << output.size() << std::endl;
    throw std::runtime_error(oss.str());
  }

}

void loadX(const std::string &filename, size_t nrow, size_t ncol, std::vector<double> &output) {
  std::ifstream X;
  X.open(filename);
  if (X.is_open()) {
    std::string line;

    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    boost::char_separator<char> sep{","};

    size_t n = 0;
    while (std::getline(X, line)) {
      ++n;
      tokenizer tok{line, sep};
      size_t m = 0;

      for (const auto &t : tok) {
        ++m;
        double temp = boost::lexical_cast<double>(t);
        output.push_back(temp);
      }

      if (m != ncol) {
        std::ostringstream oss;
        oss << "num of columns is not correct for row " << n << ", expecting " << ncol;
        throw std::runtime_error(oss.str());
      }

    }

    X.close();

    if (n != nrow) {
      std::ostringstream oss;
      oss << "num of rows is not correct for filename " << filename << ", expecting " << nrow << ", actually " << n
          << std::endl;
      throw std::runtime_error(oss.str());
    }
  }

  if ((nrow * ncol) != output.size()) {
    throw std::runtime_error("length of X is not correct! Abort!");
  }

}

void loadBounds(const std::string &filename,
                size_t nrow,
                std::vector<double> &lowerBounds,
                std::vector<double> &upperBounds) {
  std::ifstream X;
  X.open(filename);
  if (X.is_open()) {
    std::string line;

    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    boost::char_separator<char> sep{","};

    size_t n = 0;
    while (std::getline(X, line)) {
      ++n;
      tokenizer tok{line, sep};
      size_t m = 0;

      for (const auto &t : tok) {
        double temp = boost::lexical_cast<double>(t);
        if (m == 0) {
          lowerBounds.push_back(temp);
        } else if (m == 1) {
          upperBounds.push_back(temp);
        } else {
          std::ostringstream oss;
          oss << "the num of columns in the bounds file " << filename << " is not correct, m=" << m << std::endl;
          throw std::runtime_error(oss.str());
        }
        m++;
      }
    }

    if (n != nrow) {
      std::ostringstream oss;
      oss << "num of rows is not correct for filename " << filename << ", expecting " << nrow;
      throw std::runtime_error(oss.str());
    }
  }

  X.close();

}

int main(int argc, char **argv) {

  gflags::SetUsageMessage("some usage message");
  gflags::SetVersionString("1.0.0");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_verbose) {
    std::cout << gflags::ProgramInvocationShortName() << ": ";
    std::cout << FLAGS_config << std::endl;
    gflags::ShutDownCommandLineFlags();
  }

  std::cout << "reading ini file from " << FLAGS_config << std::endl;

  INIReader reader(FLAGS_config);

  if (reader.ParseError() != 0) {
    std::cout << "Can't load 'config.ini'\n";
    return 1;
  }

  std::string workingDirectory = reader.Get("task", "working_dir", "/home/cheng/");
  std::string problemKeyName = reader.Get("task", "key", "INVALID_KEY");
  int problemId = reader.GetInteger("task", "id", -1);

  std::string XFilename = workingDirectory + reader.Get("files", "X", "X.csv");
  std::string yFilename = workingDirectory + reader.Get("files", "y", "y.csv");
  std::string wFilename = workingDirectory + reader.Get("files", "w", "w.csv");
  std::string boundFilename = workingDirectory + reader.Get("files", "bound", "bound.csv");
  std::string summaryFilename = workingDirectory + reader.Get("files", "summary", "summary.csv");

  bool applyWeight = reader.GetBoolean("hyper", "use_weight", false);
  bool applyRegu = reader.GetBoolean("hyper", "reg", false);
  double minLambda = reader.GetReal("hyper", "min_lambda", 1.0);
  double maxLambda = reader.GetReal("hyper", "max_lambda", 100);
  int32_t nSample = reader.GetInteger("hyper", "nSample", 1);
  int32_t nFeature = reader.GetInteger("hyper", "nFeature", 1);

  std::cout << "we are solving problem " << problemKeyName << ", id=" << problemId << "\nX from " << XFilename
            << "\ny from " << yFilename << "\nw from " << wFilename << "\n #feature " << nFeature << "\n #sample "
            << nSample << "\n useWeight=" << applyWeight << "\n reg = " << applyRegu
            << "\n working_dir is " << workingDirectory << std::endl;

  shape_t<1> betaShape(nFeature);

  /**************************************************************
   *  create bounds
  **************************************************************/
  std::vector<double> betaLowerBoundVector;
  std::vector<double> betaUpperBoundVector;

  loadBounds(boundFilename, nFeature, betaLowerBoundVector, betaUpperBoundVector);

  assert(betaLowerBoundVector.size() == betaUpperBoundVector.size());

  std::shared_ptr<ndarray<double, 1>> betaLowerBound = std::shared_ptr<ndarray<double, 1>>(new ndarray<double, 1>(
      betaShape,
      betaLowerBoundVector.begin(),
      betaLowerBoundVector.end()));

  std::shared_ptr<ndarray<double, 1>> betaUpperBound = std::shared_ptr<ndarray<double, 1>>(new ndarray<double, 1>(
      betaShape,
      betaUpperBoundVector.begin(),
      betaUpperBoundVector.end()));

  /**************************************************************
   *  create weight vector
  **************************************************************/
  std::vector<double> wAsVector;
  size_t wDim = applyRegu ? (nSample + nFeature) : nSample;
  shape_t<1> wShape(wDim);
  loadWeight(wFilename, nSample, wAsVector);
  if (applyRegu) {
    for (int i = 0; i < nFeature; ++i) {
      wAsVector.push_back(1.0);
    }
  }

  std::shared_ptr<ndarray<double, 1>> wVector = std::shared_ptr<ndarray<double, 1>>(new ndarray<double, 1>(
      wShape,
      wAsVector.begin(),
      wAsVector.end()));

  /**************************************************************
   *  create y
  **************************************************************/
  std::vector<double> yAsVector;
  size_t yDim = applyRegu ? (nSample + nFeature) : nSample;
  shape_t<1> yShape(yDim);
  loady(yFilename, nSample, yAsVector);
  if (applyRegu) {
    for (int i = 0; i < nFeature; ++i) {
      yAsVector.push_back(0.0);
    }
  }

  std::shared_ptr<ndarray<double, 1>> yVector = std::shared_ptr<ndarray<double, 1>>(new ndarray<double, 1>(
      yShape,
      yAsVector.begin(),
      yAsVector.end()));

  /**************************************************************
   *  create X with placeholder of appended regu matrix
  **************************************************************/
  std::vector<double> XAsVector;
  loadX(XFilename, nSample, nFeature, XAsVector);

  size_t X_n = applyRegu ? (nSample + nFeature) : nSample;
  // augment the regu terms
  if (applyRegu) {
    for (size_t i = 0; i < nFeature; ++i) {
      for (size_t k = 0; k < nFeature; ++k) {
        XAsVector.push_back(0);
      }
    }
  }

  shape_t<2> XShape(X_n, nFeature);
  std::shared_ptr<ndarray<double, 2>> XMatrix = std::shared_ptr<ndarray<double, 2>>(new ndarray<double, 2>(
      XShape,
      XAsVector.begin(),
      XAsVector.end()));

  /**************************************************************
   *  create lambda sequence
  **************************************************************/
  std::vector<double> lambdaVector;
  double step = (maxLambda - minLambda) / 100.0;
  for (size_t i = 0; i < 100; ++i) {
    lambdaVector.push_back(minLambda + static_cast<double>(i) * step);
  }

  /**************************************************************
   * start solving the problem.
   * Iterate through all the lambda in the lambda sequence
   *************************************************************/
  std::ofstream summaryFile;
  summaryFile.open(summaryFilename, std::ofstream::out);

  summaryFile << "lambda,";
  for (size_t i = 0; i < nFeature - 1; ++i) {
    summaryFile << "feature_" << i << ",";
  }
  summaryFile << "feature_" << nFeature - 1 << ",total_ad,reg" << std::endl;

  size_t lambdaId = 0;
  for (auto &lambda : lambdaVector) {
    std::cout << "working on lambda = " << lambda << std::endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    Model::t M = new Model(problemKeyName.c_str());

    auto _M = finally([&]() { M->dispose(); });

    std::ostringstream filenmaeOss;
    filenmaeOss << workingDirectory << lambda << ".log";
    TaskLogger logger(filenmaeOss.str());

    //M->setLogHandler([=](const std::string &msg) { std::cout << msg << std::flush; });
    M->setLogHandler([&logger](const std::string &msg) { return logger.log(msg); });

    Variable::t beta = M->variable("beta", nFeature, Domain::inRange(betaLowerBound, betaUpperBound));
    Variable::t tau = M->variable("tau", X_n, Domain::greaterThan(0));
    //Variable::t tau = M->variable("tau", X_n, Domain::unbounded());

    // Create constraints
    // M->constraint(beta->index(1), Domain::lessThan(10.0));

    size_t rawSizeX = nSample * nFeature;
    /*
     * the magic part, update the lambda
     */
    for (auto i = 0; i < nFeature; ++i) {
      for (auto k = 0; k < nFeature; ++k) {
        if (i == k) {
          XAsVector[rawSizeX + i * nFeature + k] = lambda;
        }
      }
    }

    shape_t<2> XShape(X_n, nFeature);
    std::shared_ptr<ndarray<double, 2>> XMatrix = std::shared_ptr<ndarray<double, 2>>(new ndarray<double, 2>(
        XShape,
        XAsVector.begin(),
        XAsVector.end()));

    auto residual = Expr::sub(yVector, Expr::mul(XMatrix, beta));
    auto residualAndTau = Expr::add(residual, tau);

    // y -X*beta > -t
    M->constraint("negativeCase", residualAndTau, Domain::greaterThan(0.0));

    // y -X*beta < t
    M->constraint("positiveCase", Expr::sub(tau, residual), Domain::greaterThan(0.0));

    // beta > lowerBound
    M->constraint("betaLowerBound", Expr::sub(beta, betaLowerBound), Domain::greaterThan(0.0));

    // beta < upperBound
    M->constraint("betaUpperBound", Expr::sub(beta, betaUpperBound), Domain::lessThan(0.0));

    // Set the objective function to minimize tau
    M->objective("obj", ObjectiveSense::Minimize, Expr::dot(tau, wVector));

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << "model creation took " << time_span.count() << " seconds." << std::endl;

    // Solve the problem
    M->solve();

    // Get the solution values
    auto sol = beta->level();
    summaryFile << lambda << ",";
    summaryFile << M->primalObjValue() << ",";
    double reg = 0.0;
    for (auto beta : *sol) {
      if (beta < 1e-6) {
        summaryFile << 0.0 << ",";
      } else {
        summaryFile << beta << ",";
      }
      reg += fabs(beta) * lambda;
    }

    summaryFile << reg;
    summaryFile << std::endl;

    // Get the solution objection value
    std::cout << "total absolute deviate " << M->primalObjValue() << ", in which " << reg
              << " is from regu terms" << std::endl;
  }

  summaryFile.close();
}
