import argparse
import sys
from typing import Optional

from ecalautoctrl import QueryOMS


class FillingSchemeNameGetter:
    """
    Get the name of a injection scheme for the run with given number.
    """

    def __init__(self):
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument(
            "--run",
            dest="run_number",
            type=int,
            help="Run number.",
            required=True,
        )

        self.args = self.parser.parse_args()

    def __call__(self) -> Optional[str]:
        """
        Get and return the injection scheme name corresponding to the provided run number.
        """
        qom = QueryOMS()
        run_number = self.args.run_number
        fill_number = qom.get_run_attribute(run_number, "fill_number")

        return self.get_filling_scheme_name_from_fill_number(fill_number=fill_number)

    @staticmethod
    def get_filling_scheme_name_from_fill_number(fill_number: int) -> Optional[str]:
        """
        Get filling scheme name of a fill with provided number.
        """
        qom = QueryOMS()
        query = qom.oms.query("fills")
        query.verbose = False
        query.filter("fill_number", fill_number)
        query.attrs(["fill_number", "injection_scheme"])
        if fill_data := query.data().json().get("data"):
            return fill_data[0].get("attributes", {}).get("injection_scheme")
        return None


if __name__ == "__main__":
    scheme_name_getter = FillingSchemeNameGetter()
    result = scheme_name_getter()
    print(result)
    sys.exit(0 if result is not None else 1)
