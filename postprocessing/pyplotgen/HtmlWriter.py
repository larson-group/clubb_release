import os

class HtmlWriter:
    '''

    '''

    def __init__(self, cases):
        '''

        '''
        self.cases = cases

    def save(self):
        '''

        :return:
        '''
        print("Generating html output file")
        output_file = open("output/results.html", "w")
        self.write_header(output_file)
        for r, dirs, files in os.walk(os.curdir+"/output"):
            for folder in dirs:
                self.write_folder(output_file, folder)
        self.write_footer(output_file)
        output_file.close()


    # def write_folder(self, file, folder):
    #     '''
    #
    #     :param folder_name:
    #     :return:
    #     '''
    #     for _, __, files in os.walk(folder):
    #         file.write(
    #             '<table>\n' +
    #             '<tbody>\n' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             '' +
    #             ''
    #         )
    #         # <table style="height: 75px;" width="357">
    #         # <tbody>
    #         # <tr>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # </tr>
    #         # <tr>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # </tr>
    #         # <tr>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # <td style="width: 111.667px;">&nbsp;</td>
    #         # </tr>
    #         # </tbody>
    #         # </table>
    #         # <p style="text-align: center;"><img src="dir" alt="" /></p>
    #         # <p style="text-align: center;">&nbsp;</p>
    #
    # def write_header(self, file):
    #     '''
    #
    #     :return:
    #     '''
    #     file.write("" + \
    #                 '<html lang="en">\n' +\
    #                 '<head>\n' + \
    #                 '<meta charset="utf-8">\n' + \
    #                 '<title>Plot Results</title>\n' + \
    #                 '<meta name="description" content="Pyplotgen output">\n' + \
    #                 '<meta name="author" content="Pyplotgen">\n' + \
    #                 '<link rel="stylesheet" href="css/styles.css?v=1.0">\n' + \
    #                 '</head>\n' + \
    #                 '<body>\n' + \
    #                   "")
    #
    #
    # def write_footer(self, file):
    #     '''
    #
    #     :return:
    #     '''
    #     file.write("</body>\n")
