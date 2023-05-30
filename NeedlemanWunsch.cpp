#include <iostream>
#include <vector>
#include <algorithm>
#include <format>

class NeedlemanWunsch {
private:
    typedef struct MatrixCell {
        int score;
        enum Trace {
            UP, LEFT, DIAGONAL
        } trace;
    } MatrixCell;
    std::vector<std::vector<MatrixCell>> scoreMatrix;

    const std::string_view subject;
    const std::string_view query;
    std::string alignedSubject;
    std::string alignedQuery;

    const int matchScore;
    const int mismatchScore;
    const int gapScore;
    int totalScore = 0;

    MatrixCell calculateCell(const size_t x, const size_t y) {
        auto u = scoreMatrix[x - 1][y].score + gapScore;
        auto l = scoreMatrix[x][y - 1].score + gapScore;
        auto d = scoreMatrix[x - 1][y - 1].score + (subject[x - 1] == query[y - 1] ? matchScore : mismatchScore);

        auto score = std::max({u, l, d});
        auto direction = score == u ? MatrixCell::UP : score == l ? MatrixCell::LEFT : MatrixCell::DIAGONAL;
        return {score, direction};
    }

    void initializeScoreMatrix() {
        scoreMatrix.resize(subject.length() + 1, std::vector<MatrixCell>(query.length() + 1));

        for (int col = 0; col < subject.length() + 1; ++col) {
            scoreMatrix[col][0].score = gapScore * col;
            scoreMatrix[col][0].trace = MatrixCell::UP;
        }

        for (int row = 0; row < query.length() + 1; ++row) {
            scoreMatrix[0][row].score = gapScore * row;
            scoreMatrix[0][row].trace = MatrixCell::LEFT;
        }
    }

    void calculateScoreMatrix() {
        for (size_t x = 1; x <= subject.length(); ++x) {
            for (size_t y = 1; y <= query.length(); ++y) {
                scoreMatrix[x][y] = calculateCell(x, y);
            }
        }
    }

    void traceBestAlignment() {
        size_t x = subject.length();
        size_t y = query.length();
        alignedSubject = subject[x] + alignedSubject;
        alignedQuery = query[y] + alignedQuery;
        totalScore = scoreMatrix[x][y].score;
        while (x || y) {
            switch (scoreMatrix[x][y].trace) {
                case MatrixCell::UP:
                    alignedSubject = subject[--x] + alignedSubject;
                    alignedQuery = '-' + alignedQuery;
                    break;
                case MatrixCell::LEFT:
                    alignedSubject = '-' + alignedSubject;
                    alignedQuery = query[--y] + alignedQuery;
                    break;
                case MatrixCell::DIAGONAL:
                    alignedSubject = subject[--x] + alignedSubject;
                    alignedQuery = query[--y] + alignedQuery;
                    break;
            }
            totalScore += scoreMatrix[x][y].score;
        }
    }

public:
    NeedlemanWunsch(
            const std::string_view subject,
            const std::string_view query,
            const int matchScore = 2,
            const int mismatchScore = -1,
            const int gapScore = -2)
            : subject(subject), query(query), matchScore(matchScore), mismatchScore(mismatchScore), gapScore(gapScore) {
        initializeScoreMatrix();
        calculateScoreMatrix();
        traceBestAlignment();
    }

    void printScoreMatrix() {
        std::cout << "          ";
        for (const auto &n: query) {
            std::cout << std::format("  {}  ", n);
        }
        std::cout << "\n";
        int n = -1;
        for (const auto &column: scoreMatrix) {
            std::cout << std::format("  {}  ", n == -1 ? ' ' : subject[n]);
            for (const auto &row: column) {
                std::cout << std::format(" {:>3} ", row.score);
            }
            std::cout << "\n";
            ++n;
        }
    }

    void printBestAlignment() {
        std::cout << alignedSubject << "\n";
        for (size_t c = 0; c < alignedSubject.length() - 1; ++c) {
            std::cout << (alignedSubject[c] == alignedQuery[c] ? '|' : ' ');
        }
        std::cout << "\n" << alignedQuery << std::endl;
    }

    [[nodiscard]] std::string getSubject() const {
        return std::string(subject);
    }

    [[nodiscard]] std::string getQuery() const {
        return std::string(query);
    }

    [[nodiscard]] int getMatchScore() const {
        return matchScore;
    }

    [[nodiscard]] int getMismatchScore() const {
        return mismatchScore;
    }

    [[nodiscard]] int getGapScore() const {
        return gapScore;
    }

    [[nodiscard]] int getTotalScore() const {
        return totalScore;
    }

    [[nodiscard]] std::vector<std::vector<MatrixCell>> getScoreMatrix() const {
        return scoreMatrix;
    }
};

int main() {
    const std::string strandA = "CACGTGATCAA";
    const std::string strandB = "AGCATCGGTTG";
    const int matchScore = 2;
    const int mismatchScore = -1;
    const int gapScore = -2;

    NeedlemanWunsch NWInstance(strandA, strandB, matchScore, mismatchScore, gapScore);

    std::cout << "STRAND #1: " << NWInstance.getSubject() << "\n";
    std::cout << "STRAND #2: " << NWInstance.getQuery() << "\n\n";

    std::cout << "SCORING SCHEME:\n"
                 "- MATCH     = " << NWInstance.getMatchScore() << "\n" <<
                 "- MISMATCH  = " << NWInstance.getMismatchScore() << "\n" <<
                 "- INDEL/GAP = " << NWInstance.getGapScore() << "\n\n";

    std::cout << "MATRIX:\n";
    NWInstance.printScoreMatrix();
    std::cout << "\nALIGNMENT:\n";
    NWInstance.printBestAlignment();
    std::cout << "\nSCORE: " << NWInstance.getTotalScore() << std::endl;

    return 0;
}