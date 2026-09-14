from pathlib import Path

from .metrics.base_calculator import MetricResult


def write_output_file(results: list[MetricResult], outfile: str) -> None:
    """
    Writes the TAP results to an output file in csv format.

    Args:
        results: the list of metric results
        outfile: the path to where the results should be written.
    """
    outstr = "Metric,Value,Flag\n"
    for res in results:
        outstr += f"{res.metric_name},{res.calculated_value:.2f},{res.flag}\n"

    # VENDOR EDIT: upstream calls ab_characterisation...utils.outputs.write_file,
    # which also handles s3:// destinations. We only ever write local paths, so
    # this is a plain local write instead of vendoring the S3 helper too.
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    Path(outfile).write_text(outstr)

    return
