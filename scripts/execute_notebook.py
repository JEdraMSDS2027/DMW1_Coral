"""Execute the report in a fresh kernel using this Python interpreter.

Run from any directory: python /path/to/project/scripts/execute_notebook.py
Temporary Jupyter state stays inside the project; original notebooks are not run.
"""

import asyncio
import argparse
import json
import os
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter


def main():
    """Run every report cell and save executed notebook and an HTML reading copy."""
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Execute a project notebook in a fresh kernel.")
    parser.add_argument("notebook",nargs="?",default="coral_cleaning_preprocessing.ipynb",
                        choices=["coral_cleaning_preprocessing.ipynb","coral_pca.ipynb",
                                 "consolidated_notebook.ipynb",
                                 "consolidated_report_integrated.ipynb"])
    args=parser.parse_args()
    if sys.platform == "win32":
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    state = root / ".jupyter"
    kernel = state / "kernels" / "coral-project"
    kernel.mkdir(parents=True, exist_ok=True)
    (state / "runtime").mkdir(exist_ok=True)
    spec = {"argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
            "display_name": "Coral preprocessing", "language": "python"}
    (kernel / "kernel.json").write_text(json.dumps(spec), encoding="utf-8")
    os.environ["JUPYTER_PATH"] = str(state) + os.pathsep + os.environ.get("JUPYTER_PATH", "")
    os.environ["JUPYTER_RUNTIME_DIR"] = str(state / "runtime")
    notebook_path = root / args.notebook
    notebook = nbformat.read(notebook_path, as_version=4)
    print(f"Executing {args.notebook} in a fresh local kernel...", flush=True)
    client = NotebookClient(notebook, timeout=180, startup_timeout=30, kernel_name="coral-project",
                            resources={"metadata": {"path": str(root)}}, allow_errors=False)
    client.execute()
    # Refresh the opening abstract from this run's computed report fragment.
    if args.notebook in {"consolidated_notebook.ipynb", "consolidated_report_integrated.ipynb"}:
        abstract = None
        for cell in notebook.cells:
            if cell.cell_type != "code":
                continue
            retained_outputs = []
            for output in cell.get("outputs", []):
                if output.get("metadata", {}).get("report_fragment") == "abstract":
                    abstract = output["data"]["text/markdown"]
                else:
                    retained_outputs.append(output)
            cell.outputs = retained_outputs
        if abstract is None:
            raise ValueError("Final notebook did not generate its abstract.")
        for cell in notebook.cells:
            if cell.metadata.get("report_fragment") == "abstract":
                cell.source = "**Abstract**\n\n" + abstract
    nbformat.validate(notebook)
    nbformat.write(notebook, notebook_path)
    exporter = HTMLExporter(exclude_input_prompt=True, exclude_output_prompt=True)
    body, _ = exporter.from_notebook_node(notebook)
    if args.notebook == "consolidated_report_integrated.ipynb":
        html = root/"outputs"/"integrated_report"/"consolidated_report_integrated.html"
    elif args.notebook == "consolidated_notebook.ipynb":
        html = root/"outputs"/"final_report"/"consolidated_notebook.html"
    elif args.notebook == "coral_pca.ipynb":
        html = root/"outputs"/"pca"/"coral_pca.html"
    else:
        html = root/"processed"/"coral_cleaning_preprocessing.html"
    html.write_text(body, encoding="utf-8")
    print(f"Executed {sum(cell.cell_type == 'code' for cell in notebook.cells)} code cells successfully.")
    print(f"Report: {notebook_path}")
    print(f"HTML: {html}")


if __name__ == "__main__":
    main()
