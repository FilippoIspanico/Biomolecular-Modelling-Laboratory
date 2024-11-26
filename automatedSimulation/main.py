from src.MDRunner import MDRunner


def main():
    runner = MDRunner()
    runner.prepare_protein()
    runner.write_conf()
    runner.run()
    runner.computeQValue()
    runner.compute_rmsf(plot=True, scale=True)


if __name__ == "__main__":
    main()
