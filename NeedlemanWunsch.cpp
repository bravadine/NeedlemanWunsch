#include <iostream>
#include <vector>
#include <algorithm>
#include <format>

class NeedlemanWunsch {
public:
    std::string alignedSubject;
    std::string alignedQuery;
    int score;

    NeedlemanWunsch(
            const std::string_view subject,
            const std::string_view query,
            const int matchScore,
            const int mismatchScore,
            const int gapScore)
            : subject(subject), query(query), matchScore(matchScore), mismatchScore(mismatchScore), gapScore(gapScore) {
        initializeScoreMatrix();
        calculateScoreMatrix();
        traceBestAlignment();
    }

    ~NeedlemanWunsch() = default;

    auto getScoreMatrix() const {
        return scoreMatrix;
    }

    void viewScoreMatrix();

    void viewBestAlignment();

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

    const int matchScore;
    const int mismatchScore;
    const int gapScore;

    void calculateCell(const size_t x, const size_t y) {
        auto u = scoreMatrix[x - 1][y].score + gapScore;
        auto l = scoreMatrix[x][y - 1].score + gapScore;
        auto d = scoreMatrix[x - 1][y - 1].score + (subject[x - 1] == query[y - 1] ? matchScore : mismatchScore);

        auto cellScore = std::max({u, l, d});
        auto cellTrace = cellScore == u ? MatrixCell::UP : cellScore == l ? MatrixCell::LEFT : MatrixCell::DIAGONAL;
        scoreMatrix[x][y] = {cellScore, cellTrace};
    }

    void initializeScoreMatrix() {
        scoreMatrix.resize(subject.length() + 1, std::vector<MatrixCell>(query.length() + 1));

        for (int col = 0; col < scoreMatrix.size(); ++col) {
            scoreMatrix[col][0] = {gapScore * col, MatrixCell::UP};
        }

        for (int row = 0; row < scoreMatrix.size(); ++row) {
            scoreMatrix[0][row] = {gapScore * row, MatrixCell::LEFT};
        }
    }

    void calculateScoreMatrix() {
        for (size_t x = 1; x <= subject.length(); ++x) {
            for (size_t y = 1; y <= query.length(); ++y) {
                calculateCell(x, y);
            }
        }
    }

    void traceBestAlignment() {
        size_t x = subject.length();
        size_t y = query.length();

        alignedSubject.reserve(x + y);
        alignedQuery.reserve(x + y);

        alignedSubject += subject[x];
        alignedQuery += query[x];
        score = scoreMatrix[x][y].score;

        while (x || y) {
            switch (scoreMatrix[x][y].trace) {
                case MatrixCell::UP:
                    alignedSubject += subject[--x];
                    alignedQuery += '-';
                    break;
                case MatrixCell::LEFT:
                    alignedSubject += '-';
                    alignedQuery += query[--y];
                    break;
                case MatrixCell::DIAGONAL:
                    alignedSubject += subject[--x];
                    alignedQuery += query[--y];
                    break;
            }
            score += scoreMatrix[x][y].score;
        }

        std::reverse(alignedSubject.begin(), alignedSubject.end());
        std::reverse(alignedQuery.begin(), alignedQuery.end());
    }
};

void NeedlemanWunsch::viewScoreMatrix() {
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

void NeedlemanWunsch::viewBestAlignment() {
    std::cout << alignedSubject << "\n";
    for (size_t c = 0; c < alignedSubject.length() - 1; ++c) {
        std::cout << (alignedSubject[c] == alignedQuery[c] ? '|' : ' ');
    }
    std::cout << "\n" << alignedQuery << std::endl;
}

int main() {
    const std::string strandA = "CACGTGATCAA";
    const std::string strandB = "AGCATCGGTTG";
    const int matchScore = 2;
    const int mismatchScore = -1;
    const int gapScore = -2;

    NeedlemanWunsch NWInstance(strandA, strandB, matchScore, mismatchScore, gapScore);

    std::cout << "STRAND #1: " << strandA << "\n";
    std::cout << "STRAND #2: " << strandB << "\n\n";

    std::cout << "SCORING SCHEME:\n"
                 "- MATCH     = " << matchScore << "\n" <<
                 "- MISMATCH  = " << mismatchScore << "\n" <<
                 "- INDEL/GAP = " << gapScore << "\n\n";

    std::cout << "MATRIX:\n";
    NWInstance.viewScoreMatrix();
    std::cout << "\nALIGNMENT:\n";
    NWInstance.viewBestAlignment();
    std::cout << "\nSCORE: " << NWInstance.score << std::endl;

    return 0;
}