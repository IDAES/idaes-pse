def logs_got_warning(records, text: str = "") -> bool:
    got_warning = False
    for record in records:
        if record.levelname == "WARNING" and text in record.msg:
            got_warning = True
            break
    return got_warning

