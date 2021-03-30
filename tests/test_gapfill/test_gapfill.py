# 3rd party
from coincidence.regressions import AdvancedFileRegressionFixture
from domdf_python_tools.paths import PathPlus

# this package
from pyms.Gapfill.Function import file2dataframe


def test_file2dataframe(tmp_pathplus: PathPlus, advanced_file_regression: AdvancedFileRegressionFixture):
	area_file = PathPlus(__file__).parent / "area.csv"
	file2dataframe(area_file).to_csv(tmp_pathplus / "area.csv", index=False, na_rep="NA")
	advanced_file_regression.check_file(tmp_pathplus / "area.csv")
