from postprocessing.pyplotgen.src.OutputHandler import writeFinalErrorLog


def test_single_process_error_log_finishes_without_running_past_the_end(tmp_path):
    temporary = tmp_path / "error_temp.log"
    final = tmp_path / "error.log"
    lines = [
        f"26-08-14 18:53:36.37{i} proc_id: 462339\tmessage {i}\n"
        for i in range(3)
    ]
    temporary.write_text("".join(lines), encoding="utf-8")

    writeFinalErrorLog(str(temporary), str(final))

    assert final.read_text(encoding="utf-8") == "".join(lines)
    assert not temporary.exists()
